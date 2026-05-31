"""
    simulate_soilwater(forcing, dates; ps, state, ETi_obs, Tsoil_obs, kstep, SolveSM_fn)

单独运行土壤水分模块（不调用完整BEPS），用于快速调试和参数率定。

# Arguments
- `forcing`    : `MetSeries`，只使用 `Prcp` 字段 [mm/hr]
- `ETi_obs`    : 各层蒸散发 [m s-1]，`nlayer × ntime` 矩阵；`nothing` 时全为0
- `Tsoil_obs`  : 各层土壤温度 [°C]，`nlayer × ntime`；`nothing` 时默认5°C（不冻结）

# Returns
DataFrame，列：`date`, `z_water`, `inf`, `θ1`…`θn`
"""
function simulate_soilwater(forcing::MetSeries, dates::AbstractVector;
  ps::ParamBEPS,
  state::StateBEPS,
  ETi_obs::Union{Nothing, AbstractMatrix} = nothing,
  Tsoil_obs::Union{Nothing, AbstractMatrix} = nothing,
  kstep::Float64 = 360.0,
  SolveSM_fn = SolveSM_BEPS)

  state = deepcopy(state)
  ntime = forcing.ntime
  kloop = round(Int, 3600.0 / kstep)
  n = Int(state.n_layer)

  θ_out = zeros(n, ntime)
  inf_out = zeros(ntime)
  zw_out = zeros(ntime)

  for i = 1:ntime
    if Tsoil_obs !== nothing
      state.Tsoil_c[1:n] .= @view Tsoil_obs[:, i]
    else
      state.Tsoil_c[1:n] .= 5.0
    end

    r_rain_g_hr = forcing.Prcp[i] / 3600.0 / 1000.0 # [mm/hr] -> [m/s]

    for _ = 1:kloop
      if ETi_obs !== nothing
        state.ETi[1:n] .= @view ETi_obs[:, i]
      else
        state.ETi[1:n] .= 0.0
      end

      state.r_rain_g = r_rain_g_hr
      UpdateSoilMoisture(state, ps, kstep; fix_sm=false, SolveSM_fn)
    end

    θ_out[:, i] .= state.θ[1:n]
    inf_out[i] = state.inf
    zw_out[i] = state.z_water
  end

  df = DataFrame(:date => dates, :z_water => zw_out, :inf => inf_out)
  for j in 1:n
    df[!, Symbol("θ$j")] = θ_out[j, :]
  end

  df
end

export simulate_soilwater
