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
  if max_Fb > 3.6 # 864 mm/day, [cm h-1]
    Δt = 1.0
  elseif max_Fb > 0.36 # 86.4 mm/day
    Δt = 30.0 # seconds
  else
    Δt = 360.0
  end
  Δt
end



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
  inf_wa = max((z_water / kstep + r_rain_g), 0) * 360000.0 # [m s-1] -> [cm h-1]
  clamp(inf_wa, 0, inf_max)
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
