function UpdateSoilMoisture(soil::Soil, kstep::Float64)
  inf, inf_max = 0.0, 0.0
  Δt, total_t, max_Fb = 0.0, 0.0, 0.0

  n = soil.n_layer
  @unpack dz, f_water, K_sat, Kavg, Kmid, b,
  ψ_sat, ψ,
  θ_sat, θ, θ_prev = soil
  θ_prev .= θ # assign current soil moisture to prev

  # Max infiltration calculation
  # K_sat * (1 + (θ_sat - θ_prev) / dz * ψ_sat / θ_sat * b )
  K_sat_ms_1 = K_sat[1] / 360000.0
  inf_max = K_sat_ms_1 * (1 + (θ_sat[1] - θ_prev[1]) / dz[1] * ψ_sat[1] * b[1] / θ_sat[1])
  inf = max((soil.z_water / kstep + soil.r_rain_g), 0)
  inf = clamp(inf, 0, inf_max)

  # Ponded water after runoff. This one is related to runoff. LHe.
  # soil.z_water = (soil.z_water / kstep + soil.r_rain_g - inf) * kstep * soil.r_drainage

  @inbounds while total_t < kstep
    # 为了解决相互依赖的关系，循环寻找稳态
    # the unsaturated soil water retention. LHe
    # Hydraulic conductivity: Bonan, Table 8.2, Campbell 1974, K = K_sat*(θ/θ_sat)^(2b+3)
    for i in 1:n
      ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
      Kmid[i] = cal_K(θ[i], θ_sat[i], K_sat[i] / 360000.0, b[i]) # Hydraulic conductivity, [m/s]
    end

    # Fb, flow speed. Dancy's law. LHE.
    # check the r_waterflow further. LHE
    for i in 1:n-1
      # 不同层土壤深度不同，能否这样写？
      # K * ψ * b / (b + 3): ?
      # the unsaturated hydraulic conductivity of soil layer
      Kavg[i] = (Kmid[i] * ψ[i] + Kmid[i+1] * ψ[i+1]) / (ψ[i] + ψ[i+1]) * (b[i] + b[i+1]) / (b[i] + b[i+1] + 6) # 计算平均的一种方案？
      Q = Kavg[i] * (2 * (ψ[i+1] - ψ[i]) / (dz[i] + dz[i+1]) + 1) # z direction
      Q_max = (θ_sat[i+1] - θ[i+1]) * dz[i+1] / kstep + soil.ETi[i+1]
      Q = min(Q, Q_max)

      soil.r_waterflow[i] = Q
      max_Fb = max(max_Fb, abs(Q))
    end
    # p.r_waterflow[n] = 0

    Δt = guess_step(max_Fb) # this_step
    total_t += Δt
    total_t > kstep && (Δt -= (total_t - kstep))

    # from there: kstep is replaced by this_step. LHE
    for i in 1:n
      if i == 1
        θ[i] += (inf - soil.r_waterflow[i] - soil.ETi[i]) * Δt / dz[i]
      else
        θ[i] += (soil.r_waterflow[i-1] - soil.r_waterflow[i] - soil.ETi[i]) * Δt / dz[i]
      end
      θ[i] = clamp(θ[i], soil.θ_res[i], θ_sat[i])
    end
  end

end


# Campbell 1974, Bonan 2019 Table 8.2
@fastmath function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ = ψ_sat * (θ / θ_sat)^(-b)
  max(ψ, ψ_sat)
end

@fastmath cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:Real} =
  K_sat * (θ / θ_sat)^(2 * b + 3)

"""
[m s-1] -> 1000*[mm s-1] -> 1000*[kg m-2 s-1]
"""
# 如果流速过快，则减小时间步长
function guess_step(max_Fb)
  # this constraint is too large
  if max_Fb > 1.0e-5 # 864 mm/day
    Δt = 1.0
  elseif max_Fb > 1.0e-6 # 86.4 mm/day
    Δt = 30.0 # seconds
  else
    Δt = 360.0
  end
  Δt
end
