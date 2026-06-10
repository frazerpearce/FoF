#!/usr/bin/env python3
"""Benchmark a Python linked-cell FoF over a nested cascade of IC files.

This is deliberately not the SciPy/cKDTree implementation.  It uses a sparse
periodic linked-cell grid implemented with Python dictionaries and lists.

I/O is excluded from the timed section: each IC file is loaded before timing,
and only a timing summary is written after all trials have completed.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import gc
import math
import re
import time

import numpy as np

BOX_SIZE = 1.0
B = 0.2
MIN_GROUP_SIZE = 10
DEFAULT_INPUT_DIR = Path("data")
DEFAULT_OUTPUT_PATH = Path("outdir") / "python_linkedcell_fof_timings.txt"
DEFAULT_TRIALS = 3
IC_PATTERN = "ics_*.txt"
IC_N_RE = re.compile(r"^ics_(\d+)\.txt$")


@dataclass(slots=True)
class UnionFind:
    parent: list[int]
    rank: list[int]

    @classmethod
    def create(cls, size: int) -> "UnionFind":
        return cls(parent=list(range(size)), rank=[0] * size)

    def find(self, item: int) -> int:
        parent = self.parent
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        parent = self.parent
        rank = self.rank

        left_root = left
        while parent[left_root] != left_root:
            parent[left_root] = parent[parent[left_root]]
            left_root = parent[left_root]

        right_root = right
        while parent[right_root] != right_root:
            parent[right_root] = parent[parent[right_root]]
            right_root = parent[right_root]

        if left_root == right_root:
            return

        if rank[left_root] < rank[right_root]:
            parent[left_root] = right_root
        elif rank[left_root] > rank[right_root]:
            parent[right_root] = left_root
        else:
            parent[right_root] = left_root
            rank[left_root] += 1


def particle_count_from_name(path: Path) -> int:
    match = IC_N_RE.match(path.name)
    if match is None:
        raise ValueError(f"cannot extract particle count from {path}")
    return int(match.group(1))


def cascade_paths(input_dir: Path) -> list[Path]:
    paths = sorted(input_dir.glob(IC_PATTERN), key=particle_count_from_name)
    if not paths:
        raise FileNotFoundError(f"found no {IC_PATTERN!r} files in {input_dir}")
    return paths


def load_positions(path: Path) -> np.ndarray:
    positions = np.loadtxt(path, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(f"{path}: expected three coordinate columns")
    return positions


def cell_key(ix: int, iy: int, iz: int, n_cells: int) -> int:
    return (ix * n_cells + iy) * n_cells + iz


def decode_cell_key(key: int, n_cells: int) -> tuple[int, int, int]:
    iz = key % n_cells
    tmp = key // n_cells
    iy = tmp % n_cells
    ix = tmp // n_cells
    return ix, iy, iz


def build_sparse_cells(positions: np.ndarray, n_cells: int) -> dict[int, list[int]]:
    scaled = np.floor(positions * n_cells).astype(np.int64)
    np.clip(scaled, 0, n_cells - 1, out=scaled)

    cells: dict[int, list[int]] = {}
    for index, (ix, iy, iz) in enumerate(scaled):
        key = cell_key(int(ix), int(iy), int(iz), n_cells)
        cells.setdefault(key, []).append(index)

    return cells


def periodic_distance2(pi: np.ndarray, pj: np.ndarray) -> float:
    dx = abs(float(pi[0] - pj[0]))
    dy = abs(float(pi[1] - pj[1]))
    dz = abs(float(pi[2] - pj[2]))

    if dx > 0.5 * BOX_SIZE:
        dx = BOX_SIZE - dx
    if dy > 0.5 * BOX_SIZE:
        dy = BOX_SIZE - dy
    if dz > 0.5 * BOX_SIZE:
        dz = BOX_SIZE - dz

    return dx * dx + dy * dy + dz * dz


def find_group_counts_linked_cells(positions: np.ndarray) -> tuple[int, int]:
    n_particles = int(positions.shape[0])
    if n_particles == 0:
        return 0, 0

    linking_length = B * (BOX_SIZE**3 / n_particles) ** (1.0 / 3.0)
    linking_length2 = linking_length * linking_length

    # floor gives cell_size >= linking_length, so 27 neighbouring cells suffice.
    n_cells = max(1, int(math.floor(BOX_SIZE / linking_length)))
    cells = build_sparse_cells(positions, n_cells)

    union_find = UnionFind.create(n_particles)
    union = union_find.union
    pos = positions

    offsets = (-1, 0, 1)

    for key, members in cells.items():
        ix, iy, iz = decode_cell_key(key, n_cells)

        for dx in offsets:
            nx = (ix + dx) % n_cells
            for dy in offsets:
                ny = (iy + dy) % n_cells
                for dz in offsets:
                    nz = (iz + dz) % n_cells
                    neighbour_key = cell_key(nx, ny, nz, n_cells)

                    if neighbour_key < key:
                        continue

                    neighbour_members = cells.get(neighbour_key)
                    if neighbour_members is None:
                        continue

                    if neighbour_key == key:
                        count = len(members)
                        for left_pos in range(count - 1):
                            left = members[left_pos]
                            pi = pos[left]
                            for right in members[left_pos + 1 :]:
                                if periodic_distance2(pi, pos[right]) <= linking_length2:
                                    union(left, right)
                    else:
                        for left in members:
                            pi = pos[left]
                            for right in neighbour_members:
                                if periodic_distance2(pi, pos[right]) <= linking_length2:
                                    union(left, right)

    find = union_find.find
    group_sizes: dict[int, int] = {}
    for index in range(n_particles):
        root = find(index)
        group_sizes[root] = group_sizes.get(root, 0) + 1

    n_raw_groups = len(group_sizes)
    n_kept_groups = sum(1 for size in group_sizes.values() if size >= MIN_GROUP_SIZE)
    return n_raw_groups, n_kept_groups


def timed_find_groups(positions: np.ndarray, trials: int) -> list[tuple[int, float, int, int]]:
    results = []

    for trial in range(1, trials + 1):
        was_enabled = gc.isenabled()
        gc.disable()
        start = time.perf_counter()
        n_raw, n_kept = find_group_counts_linked_cells(positions)
        elapsed = time.perf_counter() - start
        if was_enabled:
            gc.enable()

        results.append((trial, elapsed, n_raw, n_kept))

    return results


def write_timings(
    output_path: Path,
    rows: list[tuple[str, int, int, float, int, int]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="ascii") as output:
        output.write("# input_file n_particles trial analysis_seconds n_raw_groups n_kept_groups\n")
        for input_file, n_particles, trial, elapsed, n_raw, n_kept in rows:
            output.write(
                f"{input_file} {n_particles} {trial} "
                f"{elapsed:.9e} {n_raw} {n_kept}\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"directory containing {IC_PATTERN} files (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"timing table output path (default: {DEFAULT_OUTPUT_PATH})",
    )
    parser.add_argument(
        "--trials",
        type=int,
        default=DEFAULT_TRIALS,
        help=f"number of analysis trials per IC file (default: {DEFAULT_TRIALS})",
    )
    args = parser.parse_args()

    if args.trials < 1:
        parser.error("--trials must be at least 1")

    rows: list[tuple[str, int, int, float, int, int]] = []

    for path in cascade_paths(args.input_dir):
        n_from_name = particle_count_from_name(path)
        positions = load_positions(path)
        n_particles = int(positions.shape[0])
        if n_particles != n_from_name:
            raise ValueError(
                f"{path}: filename says {n_from_name} particles, "
                f"but file contains {n_particles} rows"
            )

        print(f"Analysing {path} ({n_particles:,} particles)")
        for trial, elapsed, n_raw, n_kept in timed_find_groups(positions, args.trials):
            rows.append((path.name, n_particles, trial, elapsed, n_raw, n_kept))
            print(
                f"  trial {trial}: {elapsed:.6f} s, "
                f"{n_kept} kept groups"
            )

    write_timings(args.output, rows)
    print(f"Wrote timings to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
