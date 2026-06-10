#!/usr/bin/env python3
"""Generate nested clustered particle coordinates in a periodic 3D unit cube.

The largest catalogue is generated once. Smaller catalogues are prefixes of that
same accepted-particle stream, so ics_0004096.txt is exactly the first 4096 rows
of ics_1048576.txt, and so on.
"""

import math
from pathlib import Path
import random

MIN_POWER = 12
MAX_POWER = 20
PARTICLE_COUNTS = tuple(2**power for power in range(MIN_POWER, MAX_POWER + 1))
MAX_PARTICLES = PARTICLE_COUNTS[-1]

SEED = 12345
OUTPUT_DIR = Path("data")
LEGACY_OUTPUT_PATH = OUTPUT_DIR / "ics.txt"

FIELD_AMPLITUDE = 1.25
PEAK_TARGET = (0.0, 0.0, 0.0)
LONG_MODE_PHASE_SCATTER = 0.15
SHORT_MODE_PHASE_SCATTER = 0.30

MODES = (
    ((1, 0, 0), 1.00),
    ((0, 1, 0), 1.00),
    ((0, 0, 1), 1.00),
    ((1, 1, 0), 0.65),
    ((1, 0, 1), 0.65),
    ((0, 1, 1), 0.65),
    ((1, -1, 0), 0.65),
    ((1, 0, -1), 0.65),
    ((0, 1, -1), 0.65),
    ((1, 1, 1), 0.45),
)


def field_value(x: float, y: float, z: float, phases: list[float]) -> float:
    field = 0.0
    for ((kx, ky, kz), weight), phase in zip(MODES, phases):
        angle = 2.0 * math.pi * (kx * x + ky * y + kz * z) + phase
        field += weight * math.cos(angle)
    return field


def density(
    x: float,
    y: float,
    z: float,
    phases: list[float],
    normalization: float,
) -> float:
    return math.exp(FIELD_AMPLITUDE * field_value(x, y, z, phases) / normalization)


def generate_phases(rng: random.Random) -> list[float]:
    phases = []
    for (kx, ky, kz), _ in MODES:
        target_phase = -2.0 * math.pi * (
            kx * PEAK_TARGET[0] + ky * PEAK_TARGET[1] + kz * PEAK_TARGET[2]
        )
        wave_number_squared = kx * kx + ky * ky + kz * kz
        scatter = (
            LONG_MODE_PHASE_SCATTER
            if wave_number_squared == 1
            else SHORT_MODE_PHASE_SCATTER
        )
        phases.append(rng.gauss(target_phase, scatter) % (2.0 * math.pi))
    return phases


def generate_positions(
    n_particles: int,
    rng: random.Random,
    phases: list[float],
    normalization: float,
    maximum_density: float,
) -> list[tuple[float, float, float]]:
    positions: list[tuple[float, float, float]] = []
    append = positions.append

    while len(positions) < n_particles:
        x = rng.random()
        y = rng.random()
        z = rng.random()
        if rng.random() * maximum_density <= density(x, y, z, phases, normalization):
            append((x, y, z))

    return positions


def output_path(n_particles: int) -> Path:
    return OUTPUT_DIR / f"ics_{n_particles:07d}.txt"


def write_positions(path: Path, positions: list[tuple[float, float, float]]) -> None:
    with path.open("w", encoding="ascii") as output:
        output.writelines(f"{x:.16e} {y:.16e} {z:.16e}\n" for x, y, z in positions)


def main() -> None:
    rng = random.Random(SEED)
    phases = generate_phases(rng)

    normalization = math.sqrt(0.5 * sum(weight**2 for _, weight in MODES))
    maximum_field = sum(abs(weight) for _, weight in MODES)
    maximum_density = math.exp(FIELD_AMPLITUDE * maximum_field / normalization)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Generating {MAX_PARTICLES:,} accepted particles")
    positions = generate_positions(
        MAX_PARTICLES,
        rng,
        phases,
        normalization,
        maximum_density,
    )

    for n_particles in PARTICLE_COUNTS:
        path = output_path(n_particles)
        write_positions(path, positions[:n_particles])
        print(f"Wrote {path} ({n_particles:,} particles)")

    write_positions(LEGACY_OUTPUT_PATH, positions[:PARTICLE_COUNTS[0]])
    print(f"Wrote {LEGACY_OUTPUT_PATH} as legacy alias for {PARTICLE_COUNTS[0]:,} particles")


if __name__ == "__main__":
    main()
