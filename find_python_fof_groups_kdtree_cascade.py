#!/usr/bin/env python3
"""Benchmark Python FoF over a nested cascade of IC files.

I/O is deliberately excluded from the timed section: each IC file is loaded
before the timer starts, and the script writes only a timing summary after all
analysis trials have completed.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import gc
import re
import time

import numpy as np
from scipy.spatial import cKDTree

BOX_SIZE = 1.0
B = 0.2
MIN_GROUP_SIZE = 10
DEFAULT_INPUT_DIR = Path("data")
DEFAULT_OUTPUT_PATH = Path("outdir") / "python_fof_timings.txt"
DEFAULT_TRIALS = 3
IC_PATTERN = "ics_*.txt"
IC_N_RE = re.compile(r"^ics_(\d+)\.txt$")


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)

        if left_root == right_root:
            return

        if self.rank[left_root] < self.rank[right_root]:
            self.parent[left_root] = right_root
        elif self.rank[left_root] > self.rank[right_root]:
            self.parent[right_root] = left_root
        else:
            self.parent[right_root] = left_root
            self.rank[left_root] += 1


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


def load_positions(path: Path) -> list[tuple[float, float, float]]:
    positions = []
    with path.open("r", encoding="ascii") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected 3 columns")
            positions.append(tuple(float(field) for field in fields))
    return positions


def find_groups(positions: list[tuple[float, float, float]]) -> list[list[int]]:
    n_particles = len(positions)
    if n_particles == 0:
        return []

    linking_length = B * (BOX_SIZE**3 / n_particles) ** (1.0 / 3.0)
    position_array = np.asarray(positions)
    tree = cKDTree(position_array, boxsize=BOX_SIZE)

    union_find = UnionFind(n_particles)
    union = union_find.union
    for left, right in tree.query_pairs(linking_length, output_type="ndarray"):
        union(int(left), int(right))

    grouped: dict[int, list[int]] = defaultdict(list)
    find = union_find.find
    for index in range(n_particles):
        grouped[find(index)].append(index + 1)

    return sorted(grouped.values(), key=lambda group: (group[0], len(group)))


def timed_find_groups(
    positions: list[tuple[float, float, float]],
    trials: int,
) -> list[tuple[int, float, int, int]]:
    results = []

    for trial in range(1, trials + 1):
        was_enabled = gc.isenabled()
        gc.disable()
        start = time.perf_counter()
        groups = find_groups(positions)
        elapsed = time.perf_counter() - start
        if was_enabled:
            gc.enable()

        kept_groups = sum(1 for group in groups if len(group) >= MIN_GROUP_SIZE)
        results.append((trial, elapsed, len(groups), kept_groups))

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


def main() -> None:
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
        n_particles = len(positions)
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


if __name__ == "__main__":
    main()
