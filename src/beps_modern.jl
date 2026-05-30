"""
    simulate(forcing, lai, dates; ps, state, ...)

Run BEPS simulation using the modern `ParamBEPS + StateBEPS` API.

# Arguments
- `forcing`     : `MetSeries` with hourly meteorological data
- `lai`         : daily LAI vector (length = number of days)
- `dates`       : `DateTime` vector matching `forcing` length
- `ps`          : `ParamBEPS` model parameters
- `state`       : `StateBEPS` initial state (copied internally, external value unchanged)
- `lon`, `lat`  : longitude and latitude [°]
- `clumping`    : canopy clumping index (default 0.85)
- `SM_obs`      : observed soil moisture [m³/m³], `nlayer × ntime` matrix; when provided,
                  `UpdateSoilMoisture` is skipped and each hourly θ is prescribed from data
- `TS_obs`   : observed soil temperature [°C], `nlayer × ntime` matrix; when provided,
                  `UpdateHeatFlux` skips Tsoil update, ice_ratio still updated from obs temps
- `VARS_STATE` : state variables to save (default `DEFAULT_STATE_EXPORT`)
- `VARS_CACHE` : `LeafCache` variables to save (default `DEFAULT_CACHE_EXPORT`)

# Returns
`(df_flux, df_ET, states, caches)` — hourly flux DataFrames, a `StateSeries`,
and a `CacheSeries`
"""
function simulate(forcing::MetSeries, lai::Vector, dates::AbstractVector;
  ps::ParamBEPS, state::StateBEPS,
  lon::FT=120.0, lat::FT=20.0,
  kstep::Float64=360.0,
  fix_snowpack=true, fix_annual_Ta=true,
  SolveSM_fn=SolveSM_BEPS,
  SM_obs::Union{Nothing, AbstractMatrix}=nothing,
  TS_obs::Union{Nothing, AbstractMatrix}=nothing,
  VARS_STATE::Vector{Symbol}=DEFAULT_STATE_EXPORT,
  VARS_CACHE::Vector{Symbol}=DEFAULT_CACHE_EXPORT, kw...) where {FT<:AbstractFloat}

  state = deepcopy(state) # 避免外部状态被修改
  met = Met()
  mid_flux = Flux()
  mid_ET = ETFlux()
  cache = LeafCache()

  ntime = forcing.ntime
  fluxes = FluxSeries(; ntime)
  fluxes_ET = ETSeries(; ntime)

  SF, VF = split_vars(VARS_STATE)
  states = StateSeries(SF, VF, layer, ntime)

  CF = split_cache_vars(VARS_CACHE)
  caches = CacheSeries(CF, 4, ntime)

  Ta_annual = mean(forcing.Tair)
  jdays = dayofyear.(dates)
  hours = hour.(dates)

  fix_sm = SM_obs !== nothing
  fix_Tsoil = TS_obs !== nothing

  clumping = ps.veg.Ω
  for i = 1:ntime
    jday = jdays[i]
    hour = hours[i]

    if fix_sm
      state.θ_prev .= state.θ
      state.θ[1:5] .= @view SM_obs[:, i]
    end
    if fix_Tsoil
      state.Tsoil_p .= state.Tsoil_c
      state.Tsoil_c[1:5] .= @view TS_obs[:, i]
    end

    fill_met!(met, forcing, i) # 驱动数据
    k = ceil(Int, i / 24)
    _lai = lai[k]

    inter_prg_jl(jday, hour, lon, lat, _lai, clumping,
      met, ps, state, mid_flux, mid_ET, cache; kstep,
      fix_snowpack, fix_annual_Ta, Ta_annual, fix_sm, fix_Tsoil, SolveSM_fn)

    fluxes[i] = mid_flux
    fluxes_ET[i] = mid_ET
    save_state!(states, state, i, SF, VF)
    save_cache!(caches, cache, i, CF)
  end
  DataFrame(fluxes), DataFrame(fluxes_ET), states, caches
end

beps_modern(args...; kwargs...) = simulate(args...; kwargs...)


export simulate, beps_modern
