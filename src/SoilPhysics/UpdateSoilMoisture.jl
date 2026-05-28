# soil moisture transport is governed by infiltration, runoff, gradient diffusion, 
# gravity, and root
# extraction through canopy transpiration.  The net water applied to the
# surface layer is the snowmelt plus precipitation plus the throughfall
# of canopy dew minus surface runoff and evaporation.
# CLM3.5 uses a zero-flow bottom boundary condition.

function update_SoilWaterFrac!(f_water::Vector{Float64}, Tsoil_c::Vector{Float64}; n::Integer=5)
  # 注意：Tsoil_c 长度通常是 n，但这里循环到 n+1，需确认 Tsoil_c 实际分配长度。
  # 假设 Tsoil_c 长度足够，或者边界处理
  @inbounds for i in 1:n+1
    val_T = i <= length(Tsoil_c) ? Tsoil_c[i] : Tsoil_c[end] # 边界保护
    if val_T > 0.0
      f_water[i] = 1.0
    elseif val_T < -1.0
      f_water[i] = 0.1
    else
      f_water[i] = 0.1 + 0.9 * (val_T + 1.0)
    end
  end
end


function cal_infiltration(θ::V, dz::V,
  K_sat::V, θ_sat::V, ψ_sat::V, b::V,
  frac_water::T, z_water::T, r_rain_g::T, kstep::T) where {T<:AbstractFloat,V<:AbstractVector{T}}

  # pdv(psi, z) = [(psi_0 + z_0) - (psi_1 + z_1)] / [z_0 - z_1] = 1 + (psi_0 - psi_1) / Δz
  dψ_cm = -(θ_sat[1] - θ[1]) * ψ_sat[1] * b[1] / θ_sat[1]
  dz_cm = dz[1] * 100.0
  inf_max = frac_water * K_sat[1] * (1 + dψ_cm / dz_cm)               # [cm h-1]
  inf_wa = max(frac_water * (z_water / kstep + r_rain_g), 0) * 360000.0 # [m s-1] -> [cm h-1]
  clamp(inf_wa, 0, inf_max)
end


"""
    update_surface_water!(st, ps, kstep) -> inf

更新地表积水 `state.z_water`：处理降雨入渗和地表径流，返回入渗率 [cm/h]。
与 `UpdateSoilMoisture` 共享相同物理，但不更新 θ，适用于观测土壤水模式。
"""
function update_surface_water!(st::S, ps::P, kstep::Float64) where {
  S<:Union{StateBEPS,Soil},P<:Union{ParamBEPS,Soil}}

  n = st.n_layer
  (; θ_sat, K_sat, b, ψ_sat) = get_hydraulic(ps)
  (; θ, f_water, Tsoil_c, dz, z_water, r_rain_g) = st
  r_drainage = ps.r_drainage

  update_SoilWaterFrac!(f_water, Tsoil_c; n)
  inf = cal_infiltration(θ, dz, K_sat, θ_sat, ψ_sat, b, f_water[1], z_water, r_rain_g, kstep)
  inf_ms = inf / 360000.0 # [cm h-1] -> [m s-1]
  st.z_water = (z_water / kstep + r_rain_g - inf_ms) * kstep * r_drainage # Ponded water after runoff
  inf
end


function solve_SM_beps(st::S, ps::P, inf::Float64, kstep::Float64) where {
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
      _Δz_cm = (dz[i] + dz[i+1])/2 * 100.0
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


# 旧版本：兼容 Soil 结构体
UpdateSoilMoisture(soil::Soil, kstep::Float64) = UpdateSoilMoisture(soil, soil, kstep)


# 新版本：JAX 风格 (st, ps) 签名
function UpdateSoilMoisture(st::S, ps::P, kstep::Float64; fix_sm::Bool=false) where {
  S<:Union{StateBEPS,Soil},P<:Union{ParamBEPS,Soil}}

  n = st.n_layer
  (; θ, θ_prev, ice_ratio) = st

  θ_prev .= θ
  inf = update_surface_water!(st, ps, kstep)
  fix_sm && return # 如果 fix_sm=true，则只更新地表积水，不改变土壤水分状态

  solve_SM_beps(st, ps, inf, kstep) # 求解土壤水

  # update ice ratio
  for i in 1:n
    ice_ratio[i] *= θ_prev[i] / θ[i]
    ice_ratio[i] = min(1.0, ice_ratio[i])
  end
end


# Campbell 1974, Bonan 2019 Table 8.2
@fastmath function cal_ψ(θ::T, θ_sat::T, ψ_sat::T, b::T) where {T<:Real}
  ψ_cm = ψ_sat * (θ / θ_sat)^(-b)
  min(ψ_cm, ψ_sat) # both inputs are in [cm] and negative
end

@fastmath cal_K(θ::T, θ_sat::T, K_sat::T, b::T) where {T<:Real} =
  K_sat * (θ / θ_sat)^(2 * b + 3)

"""
[m s-1] -> 1000*[mm s-1] -> 1000*[kg m-2 s-1]
"""
# 如果流速过快，则减小时间步长
function guess_step(max_Fb)
  # this constraint is too large
  if max_Fb > 3.6 # 864 mm/day
    Δt = 1.0
  elseif max_Fb > 0.36 # 86.4 mm/day
    Δt = 30.0 # seconds
  else
    Δt = 360.0
  end
  Δt
end

# Function to calculate soil water uptake from a layer
"""
    Root Water Uptake

- `土壤蒸发`：仅发生在表层
- `植被蒸腾`：根据根系分布，耗水可能来自于土壤的每一层
"""
function Root_Water_Uptake(st::S, Trans_o::Float64, Trans_u::Float64, Evap_soil::Float64) where {
  S<:Union{StateBEPS,Soil}}

  Trans = Trans_o + Trans_u
  st.ETi[1] = Trans / ρ_w * st.w_norm[1] + Evap_soil / ρ_w
  for i in 2:st.n_layer
    st.ETi[i] = Trans / ρ_w * st.w_norm[i]
  end
end
