#!/usr/bin/env python3
"""Compare FoF group memberships without depending on group or member order."""

import argparse
from collections import defaultdict
from pathlib import Path


DEFAULT_REFERENCE_PATH = Path("outdir") / "python_members.txt"
DEFAULT_POSITIONS_PATH = Path("data") / "ics.txt"
EXPECTED_HEADER = "# group_id size member_ids..."


def particle_count(path: Path) -> int:
    count = 0
    with path.open("r", encoding="ascii") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected 3 columns")
            try:
                tuple(float(field) for field in fields)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: coordinates must be numeric"
                ) from error
            count += 1
    return count


def load_groups(path: Path, maximum_particle_id: int) -> dict[int, frozenset[int]]:
    groups = {}
    assigned_particles = {}
    with path.open("r", encoding="ascii") as input_file:
        header = input_file.readline().rstrip("\r\n")
        if header != EXPECTED_HEADER:
            raise ValueError(
                f"{path}:1: expected header {EXPECTED_HEADER!r}, found {header!r}"
            )

        for line_number, line in enumerate(input_file, start=1):
            line_number += 1
            if not line.strip():
                continue
            if line.startswith("#"):
                raise ValueError(f"{path}:{line_number}: unexpected comment line")

            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: expected group id, size, and member ids"
                )

            try:
                group_id = int(fields[0])
                expected_size = int(fields[1])
                member_ids = [int(field) for field in fields[2:]]
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: all fields must be integers"
                ) from error

            if group_id < 1:
                raise ValueError(f"{path}:{line_number}: group id must be positive")
            if group_id in groups:
                raise ValueError(f"{path}:{line_number}: duplicate group id {group_id}")
            if expected_size < 1:
                raise ValueError(f"{path}:{line_number}: group size must be positive")
            if len(member_ids) != expected_size:
                raise ValueError(
                    f"{path}:{line_number}: expected {expected_size} member ids, "
                    f"found {len(member_ids)}"
                )
            if len(set(member_ids)) != len(member_ids):
                raise ValueError(
                    f"{path}:{line_number}: a particle appears more than once in group {group_id}"
                )

            for member_id in member_ids:
                if not 1 <= member_id <= maximum_particle_id:
                    raise ValueError(
                        f"{path}:{line_number}: particle id {member_id} is outside "
                        f"the valid range 1..{maximum_particle_id}"
                    )
                if member_id in assigned_particles:
                    raise ValueError(
                        f"{path}:{line_number}: particle id {member_id} appears in "
                        f"groups {assigned_particles[member_id]} and {group_id}"
                    )
                assigned_particles[member_id] = group_id

            groups[group_id] = frozenset(member_ids)
    return groups


def groups_by_members(
    groups: dict[int, frozenset[int]],
) -> dict[frozenset[int], list[int]]:
    indexed = defaultdict(list)
    for group_id, members in groups.items():
        indexed[members].append(group_id)
    return indexed


def describe_group(group_id: int, members: frozenset[int]) -> str:
    preview = " ".join(str(member_id) for member_id in sorted(members)[:10])
    if len(members) > 10:
        preview += " ..."
    return f"group {group_id} ({len(members)} members): {preview}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("members_file", type=Path, help="members file to compare")
    parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_REFERENCE_PATH,
        help=f"reference members file (default: {DEFAULT_REFERENCE_PATH})",
    )
    parser.add_argument(
        "--positions",
        type=Path,
        default=DEFAULT_POSITIONS_PATH,
        help=f"particle positions used to determine the valid id range (default: {DEFAULT_POSITIONS_PATH})",
    )
    args = parser.parse_args()

    try:
        maximum_particle_id = particle_count(args.positions)
        reference_groups = load_groups(args.reference, maximum_particle_id)
        candidate_groups = load_groups(args.members_file, maximum_particle_id)
    except (OSError, UnicodeError, ValueError) as error:
        parser.error(str(error))

    reference_index = groups_by_members(reference_groups)
    candidate_index = groups_by_members(candidate_groups)
    missing_memberships = reference_index.keys() - candidate_index.keys()
    extra_memberships = candidate_index.keys() - reference_index.keys()

    duplicate_mismatch = {
        members
        for members in reference_index.keys() & candidate_index.keys()
        if len(reference_index[members]) != len(candidate_index[members])
    }

    if not missing_memberships and not extra_memberships and not duplicate_mismatch:
        print(
            f"Match: {len(reference_groups)} groups in {args.members_file} "
            f"have the same member sets as {args.reference}"
        )
        return 0

    print(f"Mismatch: {args.members_file} does not match {args.reference}")
    for members in sorted(missing_memberships, key=lambda item: (-len(item), min(item, default=0))):
        group_id = reference_index[members][0]
        print(f"Missing: {describe_group(group_id, members)}")
    for members in sorted(extra_memberships, key=lambda item: (-len(item), min(item, default=0))):
        group_id = candidate_index[members][0]
        print(f"Extra: {describe_group(group_id, members)}")
    for members in duplicate_mismatch:
        print(
            "Multiplicity mismatch for membership set: "
            f"reference group ids {reference_index[members]}, "
            f"candidate group ids {candidate_index[members]}"
        )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
