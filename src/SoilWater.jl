"""
    simulate_soilwater(forcing, dates; ps, state, ETi_obs, Tsoil_obs, θ1_obs, kstep, SolveSM_fn)

Run standalone soil-water simulation driven by meteorological forcing.
"""
function simulate_soilwater(forcing::MetSeries, dates::AbstractVector;
  ps::ParamBEPS,
  state::StateBEPS,
  ETi_obs::Union{Nothing, AbstractMatrix}=nothing,
  Tsoil_obs::Union{Nothing, AbstractMatrix}=nothing,
  θ1_obs::Union{Nothing, AbstractVector}=nothing,
  kstep::Float64=360.0,
  SolveSM_fn=SolveSM_BEPS)

  state = deepcopy(state)
  ntime = forcing.ntime
  n = Int(state.n_layer)
  kloop = round(Int, 3600.0 / kstep)

  ETi_obs !== nothing && size(ETi_obs, 2) != ntime &&
    error("ETi_obs dimension mismatch: expected $ntime columns, got $(size(ETi_obs, 2))")
  Tsoil_obs !== nothing && size(Tsoil_obs, 2) != ntime &&
    error("Tsoil_obs dimension mismatch: expected $ntime columns, got $(size(Tsoil_obs, 2))")
  θ1_obs !== nothing && length(θ1_obs) != ntime &&
    error("θ1_obs length mismatch: expected $ntime, got $(length(θ1_obs))")

  θ_out = zeros(Float64, ntime, n)
  z_water = zeros(Float64, ntime)

  for i in 1:ntime
    if Tsoil_obs !== nothing
      state.Tsoil_p .= state.Tsoil_c
      state.Tsoil_c[1:n] .= @view Tsoil_obs[:, i]
    end

    state.r_rain_g = forcing.Prcp[i] / 3600.0 / 1000.0 # [mm/h] -> [m/s]
    for _ in 1:kloop
      if ETi_obs !== nothing
        state.ETi[1:n] .= @view ETi_obs[:, i]
      else
        state.ETi[1:n] .= 0.0
      end
      UpdateSoilMoisture(state, ps, kstep; fix_sm=false, SolveSM_fn)
    end

    if θ1_obs !== nothing
      state.θ[1] = clamp(θ1_obs[i], ps.hydraulic.θ_res[1], ps.hydraulic.θ_sat[1])
    end

    θ_out[i, :] .= state.θ[1:n]
    z_water[i] = state.z_water
  end

  df = DataFrame(; date=dates, z_water)
  for j in 1:n
    df[!, Symbol("θ$j")] = θ_out[:, j]
  end
  df
end

export simulate_soilwater
