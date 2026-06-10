#!/usr/bin/env python3
"""Plot FoF timing scaling from available cascade timing files.

Reads timing tables such as:

    outdir/python_fof_timings.txt
    outdir/fortran_fof_timings.txt
    outdir/c_fof_timings.txt
    outdir/cpp_fof_timings.txt

Expected non-comment rows contain at least four whitespace-separated fields:

    input_file  n_particles  trial  seconds  ...

The script groups by implementation and particle number, uses the median runtime
across trials, fits

    log10(time) = intercept + slope * log10(N)

and writes a log-log scaling plot plus a summary table.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_OUTDIR = Path("outdir")
DEFAULT_PATTERN = "*_fof_timings.txt"
DEFAULT_PLOT_PATH = DEFAULT_OUTDIR / "fof_scaling.png"
DEFAULT_SUMMARY_PATH = DEFAULT_OUTDIR / "fof_scaling_fits.txt"


@dataclass(frozen=True)
class TimingRow:
    implementation: str
    n_particles: int
    trial: int
    seconds: float
    source_path: Path


@dataclass(frozen=True)
class ScalingPoint:
    implementation: str
    n_particles: int
    median_seconds: float
    mean_seconds: float
    std_seconds: float
    min_seconds: float
    max_seconds: float
    n_trials: int


@dataclass(frozen=True)
class FitResult:
    implementation: str
    slope: float
    intercept: float
    n_points: int
    n_min: int
    n_max: int


def implementation_name(path: Path) -> str:
    name = path.stem
    suffix = "_fof_timings"
    if name.endswith(suffix):
        name = name[: -len(suffix)]
    return name.replace("_", "-")


def parse_timing_file(path: Path) -> list[TimingRow]:
    rows: list[TimingRow] = []
    impl = implementation_name(path)

    with path.open("r", encoding="ascii") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least 4 columns")

            try:
                n_particles = int(fields[1])
                trial = int(fields[2])
                seconds = float(fields[3])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: expected columns 2-4 to be n, trial, seconds"
                ) from error

            rows.append(
                TimingRow(
                    implementation=impl,
                    n_particles=n_particles,
                    trial=trial,
                    seconds=seconds,
                    source_path=path,
                )
            )

    return rows


def load_timings(paths: list[Path]) -> list[TimingRow]:
    rows: list[TimingRow] = []
    for path in paths:
        rows.extend(parse_timing_file(path))
    return rows


def summarise_timings(rows: list[TimingRow]) -> list[ScalingPoint]:
    grouped: dict[tuple[str, int], list[float]] = {}
    for row in rows:
        grouped.setdefault((row.implementation, row.n_particles), []).append(row.seconds)

    points: list[ScalingPoint] = []
    for (impl, n_particles), values in grouped.items():
        array = np.asarray(values, dtype=float)
        positive = array[array > 0.0]
        if positive.size == 0:
            continue

        points.append(
            ScalingPoint(
                implementation=impl,
                n_particles=n_particles,
                median_seconds=float(np.median(positive)),
                mean_seconds=float(np.mean(positive)),
                std_seconds=float(np.std(positive, ddof=1)) if positive.size > 1 else 0.0,
                min_seconds=float(np.min(positive)),
                max_seconds=float(np.max(positive)),
                n_trials=int(positive.size),
            )
        )

    return sorted(points, key=lambda item: (item.implementation, item.n_particles))


def fit_scaling(points: list[ScalingPoint]) -> list[FitResult]:
    by_impl: dict[str, list[ScalingPoint]] = {}
    for point in points:
        by_impl.setdefault(point.implementation, []).append(point)

    fits: list[FitResult] = []
    for impl, impl_points in sorted(by_impl.items()):
        valid = [p for p in impl_points if p.n_particles > 0 and p.median_seconds > 0.0]
        if len(valid) < 2:
            continue

        log_n = np.log10([p.n_particles for p in valid])
        log_t = np.log10([p.median_seconds for p in valid])
        slope, intercept = np.polyfit(log_n, log_t, 1)

        fits.append(
            FitResult(
                implementation=impl,
                slope=float(slope),
                intercept=float(intercept),
                n_points=len(valid),
                n_min=min(p.n_particles for p in valid),
                n_max=max(p.n_particles for p in valid),
            )
        )

    return fits


def intercept_speedups_vs_python(fits: list[FitResult]) -> dict[str, float]:
    fit_by_impl = {fit.implementation: fit for fit in fits}
    python_fit = fit_by_impl.get("python")
    if python_fit is None:
        return {fit.implementation: float("nan") for fit in fits}

    return {
        fit.implementation: 10.0 ** (python_fit.intercept - fit.intercept)
        for fit in fits
    }


def write_summary(
    path: Path,
    timing_paths: list[Path],
    points: list[ScalingPoint],
    fits: list[FitResult],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fit_by_impl = {fit.implementation: fit for fit in fits}
    intercept_speedup_by_impl = intercept_speedups_vs_python(fits)

    with path.open("w", encoding="ascii") as output:
        output.write("# FoF timing scaling summary\n")
        output.write("# Timing files:\n")
        for timing_path in timing_paths:
            output.write(f"#   {timing_path}\n")
        output.write("#\n")
        output.write("# Fits use median_seconds and log10(time) = intercept + slope * log10(N)\n")
        output.write("# intercept_speedup_vs_python = 10**(python_intercept - implementation_intercept)\n")
        output.write(
            "# implementation slope intercept intercept_speedup_vs_python "
            "n_fit_points n_min n_max\n"
        )
        for fit in fits:
            speedup = intercept_speedup_by_impl.get(fit.implementation, float("nan"))
            output.write(
                f"{fit.implementation} {fit.slope:.8e} {fit.intercept:.8e} "
                f"{speedup:.8e} {fit.n_points:d} {fit.n_min:d} {fit.n_max:d}\n"
            )

        output.write("#\n")
        output.write(
            "# implementation n_particles median_seconds mean_seconds std_seconds "
            "min_seconds max_seconds n_trials fitted_slope intercept_speedup_vs_python\n"
        )
        for point in points:
            fit = fit_by_impl.get(point.implementation)
            slope_text = f"{fit.slope:.8e}" if fit is not None else "nan"
            speedup = intercept_speedup_by_impl.get(point.implementation, float("nan"))
            speedup_text = f"{speedup:.8e}" if np.isfinite(speedup) else "nan"
            output.write(
                f"{point.implementation} {point.n_particles:d} "
                f"{point.median_seconds:.8e} {point.mean_seconds:.8e} "
                f"{point.std_seconds:.8e} {point.min_seconds:.8e} "
                f"{point.max_seconds:.8e} {point.n_trials:d} {slope_text} "
                f"{speedup_text}\n"
            )


def plot_scaling(
    path: Path,
    points: list[ScalingPoint],
    fits: list[FitResult],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    by_impl: dict[str, list[ScalingPoint]] = {}
    for point in points:
        by_impl.setdefault(point.implementation, []).append(point)

    fit_by_impl = {fit.implementation: fit for fit in fits}

    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    for impl, impl_points in sorted(by_impl.items()):
        impl_points = sorted(impl_points, key=lambda item: item.n_particles)
        n = np.asarray([p.n_particles for p in impl_points], dtype=float)
        t = np.asarray([p.median_seconds for p in impl_points], dtype=float)
        t_min = np.asarray([p.min_seconds for p in impl_points], dtype=float)
        t_max = np.asarray([p.max_seconds for p in impl_points], dtype=float)

        lower = np.maximum(t - t_min, 0.0)
        upper = np.maximum(t_max - t, 0.0)

        fit = fit_by_impl.get(impl)
        label = f"{impl}"
        if fit is not None:
            label += f" (slope={fit.slope:.2f})"

        errorbar = ax.errorbar(
            n,
            t,
            yerr=np.vstack([lower, upper]),
            marker="o",
            linestyle="none",
            capsize=3,
            label=label,
        )

        point_colour = errorbar.lines[0].get_color()

        if fit is not None:
            fit_n = np.geomspace(fit.n_min, fit.n_max, 200)
            fit_t = 10.0 ** (fit.intercept + fit.slope * np.log10(fit_n))
            ax.plot(
                fit_n,
                fit_t,
                linestyle="--",
                linewidth=1.2,
                color=point_colour,
            )

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("Number of particles, N")
    ax.set_ylabel("Analysis time / s")
    ax.set_title("FoF scaling")
    ax.grid(True, which="both", linewidth=0.5, alpha=0.35)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def find_timing_files(outdir: Path, pattern: str) -> list[Path]:
    return sorted(path for path in outdir.glob(pattern) if path.is_file())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "timing_files",
        nargs="*",
        type=Path,
        help="Timing files to read. If omitted, scans --outdir for --pattern.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Directory to scan and write outputs (default: {DEFAULT_OUTDIR})",
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help=f"Glob pattern used when timing files are omitted (default: {DEFAULT_PATTERN})",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=DEFAULT_PLOT_PATH,
        help=f"Output plot path (default: {DEFAULT_PLOT_PATH})",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY_PATH,
        help=f"Output fit summary path (default: {DEFAULT_SUMMARY_PATH})",
    )
    args = parser.parse_args()

    timing_files = args.timing_files or find_timing_files(args.outdir, args.pattern)
    if not timing_files:
        parser.error(f"no timing files found in {args.outdir} matching {args.pattern!r}")

    try:
        rows = load_timings(timing_files)
        points = summarise_timings(rows)
        fits = fit_scaling(points)
    except (OSError, UnicodeError, ValueError) as error:
        parser.error(str(error))

    if not points:
        parser.error("no positive timing measurements found")

    write_summary(args.summary, timing_files, points, fits)
    plot_scaling(args.plot, points, fits)

    print(f"Read {len(rows)} timing rows from {len(timing_files)} file(s)")
    print(f"Wrote {args.summary}")
    print(f"Wrote {args.plot}")

    intercept_speedup_by_impl = intercept_speedups_vs_python(fits)
    for fit in fits:
        speedup = intercept_speedup_by_impl.get(fit.implementation, float("nan"))
        speedup_text = f", intercept speedup vs python={speedup:.3g}" if np.isfinite(speedup) else ""
        print(
            f"{fit.implementation}: slope={fit.slope:.4f} "
            f"over N={fit.n_min}..{fit.n_max} ({fit.n_points} points)"
            f"{speedup_text}"
        )

    implementations_without_fit = sorted(
        {point.implementation for point in points} - {fit.implementation for fit in fits}
    )
    for impl in implementations_without_fit:
        print(f"{impl}: insufficient positive N points for a fit", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
