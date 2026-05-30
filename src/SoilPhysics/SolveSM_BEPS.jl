function SolveSM_BEPS(st::S, ps::P, inf::Float64, kstep::Float64) where {
  S<:Union{StateBEPS,Soil},P<:Union{ParamBEPS,Soil}}

  n = st.n_layer
  (; θ_sat, K_sat, ψ_sat, b, θ_res) = get_hydraulic(ps)
  (; dz, f_water, Kavg, Kmid, ψ, θ, ETi, r_waterflow) = st

  total_t, max_Fb = 0.0, 0.0
  @inbounds while total_t < kstep
    # 为了解决相互依赖的关系，循环寻找稳态
    # the unsaturated soil water retention. LHe
    # Hydraulic conductivity: Bonan, Table 8.2, Campbell 1974, K = K_sat*(θ/θ_sat)^(2b+3)
    for i in 1:n
      ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
      Kmid[i] = f_water[i] * cal_K(θ[i], θ_sat[i], K_sat[i], b[i]) # Hydraulic conductivity, [cm/h]
    end

    # Fb, flow speed. Dancy's law. LHE.
    # check the r_waterflow further. LHE
    for i in 1:n-1
      # 不同层土壤深度不同，能否这样写？
      # K * ψ * b / (b + 3): ?
      # the unsaturated hydraulic conductivity of soil layer
      Kavg[i] = (Kmid[i] * ψ[i] + Kmid[i+1] * ψ[i+1]) / (ψ[i] + ψ[i+1]) * (b[i] + b[i+1]) / (b[i] + b[i+1] + 6) # 计算平均的一种方案？
      # [(ψ[i] + z_i) - (ψ[i+1] + z_i+1)] / (z_i - z_i+1) = 1 - (ψ[i+1] - ψ[i]) / Δz
      _Δz_cm = (dz[i] + dz[i+1]) / 2 * 100.0
      grad_ψ = 1 - (ψ[i+1] - ψ[i]) / _Δz_cm
      Q = Kavg[i] * grad_ψ # [cm h-1]

      # `Q_max`出现了单位不匹配的问题，导致Q_max未发挥作用
      Q_max = ((θ_sat[i+1] - θ[i+1]) * dz[i+1] / kstep + ETi[i+1]) * 360000.0 # [m s-1] -> [cm h-1]
      Q = min(Q, Q_max)

      r_waterflow[i] = Q
      max_Fb = max(max_Fb, abs(Q))
    end
    # p.r_waterflow[n] = 0

    Δt = guess_step(max_Fb) # this_step
    total_t += Δt
    total_t > kstep && (Δt -= (total_t - kstep))
    inf_ms = inf / 360000.0

    # from there: kstep is replaced by this_step. LHE
    for i in 1:n
      if i == 1
        θ[i] += (inf_ms - r_waterflow[i] / 360000.0 - ETi[i]) * Δt / dz[i] # [cm h-1] -> [m s-1]
      else
        θ[i] += ((r_waterflow[i-1] - r_waterflow[i]) / 360000.0 - ETi[i]) * Δt / dz[i]
      end
      θ[i] = clamp(θ[i], θ_res[i], θ_sat[i])
    end
  end
end
