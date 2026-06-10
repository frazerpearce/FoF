#!/usr/bin/env python3
"""Plot x-y positions for particles in retained FoF groups."""

from pathlib import Path

import matplotlib.pyplot as plt


POSITIONS_PATH = Path("data") / "ics.txt"
GROUPS_PATH = Path("outdir") / "python_groups.txt"
MEMBERS_PATH = Path("outdir") / "python_members.txt"
OUTPUT_PATH = Path("outdir") / "python_groups_xy.png"
PLOT_HALF_WIDTH = 0.2


def centered_coordinate(coordinate: float) -> float:
    return (coordinate + 0.5) % 1.0 - 0.5


def load_positions(path: Path) -> list[tuple[float, float, float]]:
    positions = []
    with path.open("r", encoding="ascii") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected 3 columns")
            positions.append(tuple(float(field) for field in fields))
    return positions


def load_members(path: Path) -> dict[int, list[int]]:
    groups: dict[int, list[int]] = {}
    with path.open("r", encoding="ascii") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            if line.startswith("#") or not line.strip():
                continue

            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: expected group id, size, and member ids"
                )

            group_id = int(fields[0])
            expected_size = int(fields[1])
            member_ids = [int(field) for field in fields[2:]]
            if len(member_ids) != expected_size:
                raise ValueError(
                    f"{path}:{line_number}: expected {expected_size} member ids"
                )
            groups[group_id] = member_ids
    return groups


def load_centers(path: Path) -> dict[int, tuple[float, float, float]]:
    centers = {}
    with path.open("r", encoding="ascii") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            if line.startswith("#") or not line.strip():
                continue

            fields = line.split()
            if len(fields) != 5:
                raise ValueError(
                    f"{path}:{line_number}: expected group id, size, and 3 center coordinates"
                )
            centers[int(fields[0])] = tuple(float(field) for field in fields[2:])
    return centers


def main() -> None:
    positions = load_positions(POSITIONS_PATH)
    groups = load_members(MEMBERS_PATH)
    centers = load_centers(GROUPS_PATH)

    missing_centers = sorted(set(groups) - set(centers))
    if missing_centers:
        raise ValueError(f"Missing centers for group ids: {missing_centers}")

    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    cmap = plt.get_cmap("tab20")

    for color_index, group_id in enumerate(sorted(groups)):
        member_ids = groups[group_id]
        xs = [centered_coordinate(positions[member_id - 1][0]) for member_id in member_ids]
        ys = [centered_coordinate(positions[member_id - 1][1]) for member_id in member_ids]
        ax.scatter(
            xs,
            ys,
            s=36,
            color=cmap(color_index % cmap.N),
            label=f"group {group_id} (n={len(member_ids)})",
            edgecolors="black",
            linewidths=0.35,
        )
        center_x, center_y, _ = centers[group_id]
        ax.scatter(
            centered_coordinate(center_x),
            centered_coordinate(center_y),
            s=150,
            marker="*",
            color=cmap(color_index % cmap.N),
            edgecolors="black",
            linewidths=1.0,
            zorder=4,
        )

    ax.set_xlim(-PLOT_HALF_WIDTH, PLOT_HALF_WIDTH)
    ax.set_ylim(-PLOT_HALF_WIDTH, PLOT_HALF_WIDTH)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("FoF group members and centers in x-y")
    ax.grid(True, alpha=0.25)

    if not groups:
        ax.text(0.0, 0.0, "No retained groups", ha="center", va="center")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PATH, dpi=200)
    plt.close(fig)

    print(f"Plotted {sum(len(members) for members in groups.values())} particles")
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
