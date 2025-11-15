# ES313 — Simulation & Modeling in Julia

This repository contains the modules, scripts, and tools developed for the **ES313 Simulation Project**, focusing on orbital mechanics, satellite constellation modeling, visibility analysis, and Earth-coverage computation using Julia.  
The objective is to provide a clean, modular, and reusable framework for simulating **LEO constellations**, computing coverage performance, and generating visualisations.

## Repository Structure

### 📁 Core Julia Modules
- `SatsLEO.jl` — Orbital mechanics, ECI/ECEF transformations, Walker constellation generation, visibility checks, and coverage evaluation.

### 📁 Notebooks / Scripts
- Satellite constellation generation and visualization  
- Coverage fraction computations  
- Coverage heatmaps over latitude–longitude grids  

## 📦 Dependencies
Primary Julia packages required:
- `Plots.jl`  
- `LinearAlgebra`

## 🛰️ Example Usage (LEO Constellations)
```julia
using .SatsLEO             # Import the local SatsLEO module (must be in the same project)

P, S, F = 4, 12, 1         # Walker-Delta parameters: P = number of planes, S = sats per plane, F = phasing
i_deg = 50                 # Orbital inclination (degrees)
h_km = 800                 # Altitude above Earth (km)
a = Re + h_km * 1e3        # Semi-major axis = Earth radius + altitude

sats = walker_delta(P, S, F, i_deg, a)        # Construct the Walker constellation (myconstellation also available)
mean_coverage_fraction(sats, -i_deg, i_deg, 10)   # Compute mean coverage between ±inclination with 10° elevation mask

```

## 📜 License
MIT License

## 🤝 Contributions
Contribution of OpenAI ChatGPT 5.1
