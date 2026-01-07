# OrbitalConstellations.jl

## Introduction

This project was developed in the context of the course **ES313 – Mathematical Modeling and Computer Simulation**
at the **Royal Military Academy (Belgium)**.

The objective of this module is to model, analyze, visualize, and optimize satellite constellations
(mainly LEO) with respect to Earth coverage, constellation geometry, number of satellites, and PDOP performance.

The codebase is written in **Julia**.  
**All comments and docstrings inside the source code are written in French.**

---

## Overview

The `OrbitalConstellations` module provides tools to:

- Define satellites and orbital planes  
- Build Walker and generalized Walker constellations  
- Compute satellite positions in ECI and ECEF frames  
- Evaluate instantaneous and average Earth coverage  
- Optimize constellations using genetic algorithms  
- Optimize constellations with a fixed number of satellites  
- Optimize constellations using PDOP-based criteria  
- Visualize constellations and coverage maps  

---

## Installation / Usage

Place all source files in the same directory and load the module:

```julia
include("OrbitalConstellations.jl")
using .OrbitalConstellations
```

---

## Basic Example: Create a Walker Constellation

```julia
P = 6
S = 4
F = 1
i_deg = 55.0
a = Re + 550e3

sats = walker_delta(P, S, F, i_deg, a)
```

This creates a Walker-Delta constellation with 24 satellites.

---

## Compute Coverage

```julia
eps_deg = 10.0
vec = [4 5 4 5]

cov, N = eval_constellation(vec, F, i_deg, a, eps_deg)
```

This evaluates the average Earth coverage and returns the total number of satellites.

---

## Optimization (Variable Number of Satellites)

```julia
best_vec, cov_final, N_final =
    evolve_vec(
        8,
        16,
        1,
        55.0,
        Re + 600e3,
        10.0;
        generations = 50
    )
```

This runs a genetic algorithm to maximize coverage while controlling the constellation size.

---

## Optimization (Fixed Number of Satellites)

```julia
best_vec, cov =
    evolve_vec_fixedN(
        8,
        20,
        1,
        55.0,
        Re + 600e3,
        10.0
    )
```

This optimizes the distribution of a fixed number of satellites across orbital planes.

---

## PDOP-Based Optimization

```julia
best_vec, mpdop, cov, N =
    evolve_vec_pdop(
        8,
        16,
        1,
        55.0,
        Re + 600e3,
        10.0
    )
```

This optimization prioritizes constellation geometry quality using PDOP metrics.

---

## Visualization

### 3D Constellation View

```julia
plot_constellation(sats, 0.0)
```

### Coverage Heatmap

```julia
show_coverage_heatmap(sats, 0.0, 10.0)
```

---

## Notes

- Multi-threading is used extensively for performance.  
- Caching mechanisms are implemented to accelerate genetic algorithms.  
- The Earth is assumed spherical with constant rotation.  
- Orbits are assumed circular. 
- AI assistance (OpenAI GPT-5.1 / GPT-5.2) was used for documentation and code review.

---

## Course Context

This project is purely academic and was developed as part of  
**ES313 – Mathematical Modeling and Computer Simulation**  
at the **Royal Military Academy (Belgium)**.