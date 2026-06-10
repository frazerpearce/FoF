# FoF: A Multi-Language Friends-of-Friends Benchmark

## Overview

This project implements the same Friends-of-Friends (FoF) group-finding algorithm in multiple programming languages and compares their scaling, performance, and development styles.

The original motivation was not to identify the fastest language, but to explore how scientific programming has evolved over the last forty years:

* BBC BASIC
* Pascal
* Fortran
* C
* C++
* Python

Each implementation reads identical initial conditions and produces identical group catalogues. Timing cascades are then used to compare scaling behaviour over a range of problem sizes.

The project also serves as a small historical tour through the development of scientific computing, from hand-written algorithms in early microcomputer languages to modern library-driven workflows.

---

## The Test Problem

The benchmark uses a standard Friends-of-Friends algorithm applied to synthetic particle distributions.

A cascade of increasingly large particle catalogues is generated:

The benchmark uses particle counts from 2¹² to 2²⁰ in powers of two.

The catalogues are nested. The first 2¹² particles of the largest catalogue form the 2¹² dataset, the first 2¹³ particles form the 2¹³ dataset, and so on. This ensures that every implementation is solving exactly the same problem at each particle count.

---

## Implementations

### BBC BASIC

A deliberately naïve implementation using an O(N²) all-pairs search.

This version exists primarily as a historical reference point and demonstrates the importance of algorithmic complexity.

### Pascal

A linked-cell implementation using structured programming techniques typical of the late 1980s and early 1990s.

Interestingly, Pascal performs extremely well despite largely disappearing from modern scientific computing.

### Fortran

A traditional high-performance linked-cell implementation using simple arrays and explicit loops.

This is typically the fastest implementation in the benchmark.

### C

A pointer-based linked-cell implementation.

### C++

A class-based linked-cell implementation.

### Python (linked-cell)

A direct translation of the linked-cell algorithm into pure Python.

This demonstrates the cost of performing large numbers of operations inside the Python interpreter.

### Python (cKDTree)

Uses SciPy's highly optimised `cKDTree` implementation.

Although written in Python, most of the heavy lifting is performed by compiled code within SciPy.

---

## What This Project Demonstrates

The benchmark illustrates several distinct effects.

### 1. Algorithms matter more than languages

The difference between O(N²) and O(N) dominates the comparison.

A poorly chosen algorithm in a "fast" language can easily be slower than a good algorithm in a "slow" language.

### 2. Ecosystems matter

Modern Python is successful largely because of its ecosystem.

The `cKDTree` implementation is competitive not because Python itself is fast, but because it provides convenient access to highly optimised compiled libraries.

### 3. Development styles have evolved

Scientific programming has shifted from:

> Write everything yourself.

to:

> Assemble and orchestrate high-quality existing components.

The benchmark attempts to capture that transition.

---

## Typical Results

The linked-cell and k-d tree implementations all exhibit approximately linear scaling:

Time ∝ N

The BBC BASIC implementation exhibits approximately quadratic scaling:

Time ∝ N²

This difference becomes dramatic at large particle counts.

---

## Files

### Initial condition generation

* `generate_nested_ics.py`

### Timing cascades

* `*_cascade.*`

### Analysis

* `plot_fof_scaling.py`

### Outputs

* `fof_scaling.png`

---

## Caveats

This is not intended to be a production-quality halo finder.

The benchmark is intentionally simple and exists primarily to compare:

* language styles
* implementation approaches
* algorithmic choices
* ecosystem effects

rather than to maximise absolute performance.

---

## Conclusion

The most important lesson from this project is that the largest performance gains over the last forty years have come not from programming languages themselves, but from improvements in algorithms, libraries, and reusable software ecosystems.

Or, put more bluntly:

> Most of the speedup came from learning not to do billions of unnecessary distance calculations.

## Licence

Copyright (c) 2026 Frazer Pearce

Released under the MIT Licence. See the LICENSE file for details.

## Disclaimer

This software is provided "as is", without warranty of any kind,
express or implied, including but not limited to the warranties of
merchantability, fitness for a particular purpose and noninfringement.

The authors shall not be liable for any claim, damages or other
liability arising from the use of this software.

This project is intended for educational and benchmarking purposes.
