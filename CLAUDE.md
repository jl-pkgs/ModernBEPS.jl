# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Identity

- **Package name:** `BEPS` (registered as `BEPS.jl`)
- **Purpose:** Boreal Ecosystem Productivity Simulator — coupled soil-plant-atmosphere (SPAC) model simulating carbon, water, and energy fluxes at hourly resolution.
- **Two execution paths:** pure-Julia (`inter_prg_jl`) and C-binding via `libbeps` shared library (`inter_prg_c`).

## Commands

```bash
# Run full test suite
julia --project test/runtests.jl

# Run a single test file
julia --project test/test-beps_modern.jl

# Start Julia REPL with project activated
julia --project

# Instantiate dependencies
julia --project -e "using Pkg; Pkg.instantiate()"
```

## Architecture

### Module Layout

```
src/
├── BEPS.jl                        # Entry: exports, libbeps loading, includes
├── BEPS_modules.jl                # Physics module exports + includes
├── beps_main.jl                   # Legacy API: beps_main (C/Julia dispatch)
├── beps_modern.jl                 # Modern API: simulate(forcing, lai, dates; ps, state)
├── inter_prg.jl                   # Julia hourly time-step orchestrator (inter_prg_jl)
├── aerodynamic_conductance.jl     # Aerodynamic resistance (V1)
├── aerodynamic_conductance_V2.jl  # Aerodynamic resistance (V2, updated)
├── heat_H_and_LE.jl               # Sensible + latent heat fluxes
├── evaporation_soil.jl            # Soil evaporation
├── evaporation_canopy.jl          # Canopy evaporation/interception
├── netRadiation.jl                # Net radiation partitioning
├── photosynthesis.jl              # Farquhar photosynthesis (coupled)
├── photosynthesis_helper.jl       # Vcmax/Jmax temperature response
├── rainfall_stage.jl              # Rainfall interception stages
├── snowpack.jl                    # Snow accumulation/melt stages
├── surface_temperature.jl         # Surface/leaf temperature iteration
├── Optim.jl                       # Objective functions + SCEUA wrapper
├── clang/                         # C library bindings (via Clang.jl)
│   ├── BEPS_c.jl                  # Entry: loads libbeps, exports C wrappers
│   ├── SOIL_c.jl                  # Soil_c struct (C-layout)
│   ├── snowpack_stage.jl          # C snowpack stage wrappers
│   └── struct_SOIL.jl             # Cstruct helpers
├── DataType/                      # All parameter and state types
│   ├── DataType.jl                # Entry: includes all DataType files
│   ├── BEPS_Param.jl              # ParamBEPS{FT,N,H,T}: top-level param struct
│   ├── BEPS_State.jl              # StateBEPS + SnowLand: mutable state structs
│   ├── CanopyLayer.jl             # Leaf, Radiation structs (overstory/understory)
│   ├── Constant.jl                # Physical constants
│   ├── LeafCache.jl               # LeafCache: per-leaf intermediate variables
│   ├── Met.jl                     # Met: hourly meteorological forcing struct
│   ├── OUTPUT.jl                  # Flux, ETFlux: output flux structs
│   ├── PhotoConsts.jl             # Photosynthesis constants
│   ├── AeroConsts.jl              # Aerodynamic constants
│   ├── StateSeries.jl             # StateSeries / CacheSeries: time-series containers
│   ├── macro.jl                   # @bounds macro and helpers
│   ├── setup.jl                   # setup() / setup_model(): JAX-style init
│   └── Params/
│       ├── Params.jl              # Entry: ParamVeg, imports ModelParams types
│       ├── ParamPhoto.jl          # Farquhar photosynthesis params (standalone)
│       ├── GlobalData.jl          # Soil/veg lookup tables
│       └── Param_Init.jl          # InitParam_Veg, InitParam_Soil, Init_Soil_T_θ!
├── SoilPhysics/
│   ├── SoilPhysics.jl             # Entry: includes soil physics
│   ├── UpdateSoilMoisture.jl      # Soil moisture update (delegates to SoilDiffEqs.jl)
│   ├── UpdateHeatFlux.jl          # Soil heat flux update
│   └── soil_water_factor_v2.jl    # Soil water stress factor
├── SPAC/
│   ├── SPAC.jl                    # Entry: includes SPAC helpers
│   ├── BEPS_helper.jl             # fill_met!, par2theta, split_vars
│   ├── Leaf.jl                    # Leaf temperature / stomatal conductance
│   ├── VCmax.jl                   # Vcmax from leaf N
│   ├── lai2.jl                    # LAI partitioning (overstory/understory)
│   ├── snow_density.jl            # Snow density model
│   └── ultilize.jl                # find_VegType, find_SoilType
└── standalone/
    └── Photosynthesis/
        └── photosynthesis.jl      # Standalone Farquhar model (no canopy coupling)
```

### Key Type Hierarchy

| Type                  | Role                                                                                                     |
| --------------------- | -------------------------------------------------------------------------------------------------------- |
| `ParamBEPS{FT,N,H,T}` | Top-level parameter container: `hydraulic::HydraulicProfile`, `thermal::ThermalProfile`, `veg::ParamVeg` |
| `StateBEPS`           | Mutable model state: soil θ, Tsoil, snow, root fraction, water factors                                   |
| `ParamVeg{FT}`        | Vegetation parameters: LAI, albedo, stomatal slope, Vcmax25, clumping                                    |
| `MetSeries`           | Hourly forcing: Tair, SWin, LWin, RH, wind, precip arrays                                                |
| `Met`                 | Single-timestep meteorological forcing (filled from MetSeries)                                           |
| `Flux`                | Hourly output: GPP, NPP, SH, LH, snow depth, water depth                                                 |
| `ETFlux`              | Detailed ET partitioning: Trans_o/u, Eil, EiS, Evap_soil                                                 |
| `LeafCache`           | Per-leaf intermediate cache (temperatures, conductances)                                                 |
| `StateSeries`         | Time-series container for selected state variables                                                       |

### Key Data Flow

1. **Setup:** `state, ps = setup(VegType, SoilType)` constructs `StateBEPS` + `ParamBEPS` from lookup tables. `ps.hydraulic` (Campbell retention) and `ps.thermal` delegate to `SoilDiffEqs.jl` types.
2. **Forcing:** `MetSeries` holds all hourly inputs. `fill_met!(met, forcing, i)` copies timestep `i` into the single-step `Met` struct.
3. **Per-timestep:** `inter_prg_jl(jday, hour, lon, lat, lai, Ω, met, ps, state, flux, etflux, cache)` runs one hour: radiation → aerodynamics → leaf temperature → photosynthesis → evapotranspiration → rainfall/snowpack → soil moisture/heat.
4. **Soil physics:** `UpdateSoilMoisture!` delegates to `SoilDiffEqs.jl` solvers; `UpdateHeatFlux!` handles heat diffusion.
5. **Output collection:** `simulate` returns `(df_flux, df_ET, states::StateSeries, caches::CacheSeries)`.

### `simulate` API (modern)

```julia
state, ps = setup("evergreen_needleleaf", "silty_clay_loam")
df_flux, df_ET, states, caches = simulate(forcing, lai, dates;
  ps, state,
  lon=120.0, lat=40.0,
  SM_obs=nothing,   # prescribe soil moisture (skips UpdateSoilMoisture)
  TS_obs=nothing,   # prescribe soil temperature (skips Tsoil update)
  VARS_STATE=[:θ, :Tsoil_c, :z_snow],
  VARS_CACHE=[:Gs_o, :Gs_u])
```

### Legacy `beps_main` API

```julia
df_flux, df_ET, states = beps_main(forcing, lai, dates;
  lon=120.0, lat=20.0,
  VegType="default", SoilType="silty_clay_loam",
  version="julia")   # or "c" to call libbeps
```

### C Library Integration

`libbeps` is loaded at module init from `deps/` or via `LazyArtifacts`. C functions are wrapped in `src/clang/BEPS_c.jl`. The `Soil_c` struct mirrors the C layout. Use `version="c"` in `beps_main` or `inter_prg_c` directly to invoke the C path.

## Testing Notes

- Tests are in `test/`. Run a single file with `julia --project test/test-beps_modern.jl`.
- `test-beps_modern.jl`: end-to-end `simulate` with the modern API.
- `test-beps_main.jl`: legacy `beps_main` consistency test (C vs Julia).
- `test-photosynthesis_standalone.jl`: standalone Farquhar model (no canopy coupling).
- `test-soil_sm.jl`: soil moisture integration tests.
