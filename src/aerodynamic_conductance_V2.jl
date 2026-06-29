# Psi_h
function cal_Ψh(ξ::FT, L::FT) where {FT<:AbstractFloat}
  # Bonan 2019, Eq 6.47
  if L >= 0
    Ψ_h = -5 * ξ
  else
    x::FT = (1 - 16 * ξ)^0.25
    Ψ_h = 2 * log((1 + x^2) / 2)
  end
  return Ψ_h
end

# phi_h
function cal_ϕh(ξ::FT, L::FT) where {FT<:AbstractFloat}
  # Bonan 2019, Eq 6.38
  if L > 0
    ϕ = 1 + 5 * ξ
  else
    ϕ = (1 - 16 * ξ)^(-0.5)
  end
  return ϕ
end


function cal_Nu(u::FT, nu_lower::FT)::FT
  # lw::T = 0.3 # leaf characteristic width =0.3 for BS
  # sigma::T = 5 # shelter factor =5 for BS
  # Re = (ud * lw / sigma) / nu_lower
  Re::FT = (u * 0.1) / nu_lower # Reynold's number
  Nu::FT = 1.0 * Re^0.5 # Nusselt number
  Nu
end

function ra_leaf_boundary(Tair::FT, u::FT)::FT
  nu_lower::FT = (13.3 + Tair * 0.07) / 1000000  # viscosity (cm2/s)
  alfaw::FT = (18.9 + Tair * 0.07) / 1000000

  Nu = cal_Nu(u, nu_lower)
  rb = min(40, 0.5 * 0.1 / (alfaw * Nu)) # leaf boundary resistance, [s m-1], `rb_o` or `rb_u`
  return rb
end


function aerodynamic_conductance_jl(h::FT, h_u::FT, z_wind::FT, clumping::FT,
  Tair::FT, u::FT, H::FT, lai_o::FT, lai_u::FT=0.0)

  k::FT = 0.4   # von Karman's constant
  cp::FT = 1010 # specific heat of air (J/kg/K)
  ρ_air::FT = 1.225 # density of air at 15 C (kg/m3)
  g::FT = 9.8 # gravitational acceleration (m/s2)

  # if wind_sp == 0.0
  G_o_a = 1 / 200.0
  G_o_b = 1 / 200.0
  G_u_a = 1 / 200.0
  G_u_b = 1 / 200.0
  ra_g = 300.0
  ra_u = 0.0
  ra_o = 0.0

  u == 0 && return ra_o, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b

  ## this is else
  d::FT = 0.8 * h   # displacement height (m)
  z0m::FT = 0.08 * h  # momentum roughness length; ≈0.1h (Chen 1999)
  z0h::FT = 0.1 * z0m               # heat roughness length; kB⁻¹ ≈ 2.3

  # 确保参考高度在 d+z0m 之上（再分析10m风 + 高冠层保护）
  z_eff, _ = safe_wind_ref(z_wind, h; d_frac=FT(0.8), z0m_frac=FT(0.08))

  ustar::FT = u * k / log((z_eff - d) / z0m) # TODO: Ψ_m is ignored here, friction velocity (m/s)

  L::FT = -(ρ_air * cp * (Tair + 273.3) * ustar^3) / (k * g * H)
  L = clamp(L, -1e6, 1e6)  # guard against H≈0
  ξ::FT = clamp((z_eff - d) / L, -2.0, 1.0)  # Bonan Figure 6.3

  Ψ_h = cal_Ψh(ξ, L)

  ra_o::FT = 1 / (k * ustar) * (log((z_eff - d) / z0h) - Ψ_h)
  ra_o = clamp(ra_o, 2.0, 500.0)

  # Bonan 2019, Eq 6.38
  ψ = cal_ϕh(ξ, L)
  ψ = min(10.0, ψ)

  #******** Leaf boundary layer resistance ******************/
  # Wind speed at tree top */
  u_h::FT = 1.1 * ustar / k  # wind speed at height h
  Le = lai_o * clumping

  ## 1. overstory
  gamma_o_m = (0.167 + 0.179 * u_h) * Le^(1.0 / 3.0)
  u_o = u_h * exp(-gamma_o_m * (1.0 - d / h))  # u(d) taking as the mean wind speed inside a stand */
  rb_o = ra_leaf_boundary(Tair, u_o)

  ## 2. Understory (林下层与树干空间)
  gamma_o_h = 0.1 + lai_o^0.75   # h ~ h_u之间主林冠下部及树干对湍流的阻滞 (使用 lai_o)
  gamma_u_m = 0.1 + lai_u^0.75     # h_u ~ 0之间林下层对湍流的阻滞 (使用 lai_u)

  u_hu = u_h * exp(-gamma_o_m * (1.0 - h_u / h)) # 风速衰减至林下层顶部 (h_u)
  u_u = u_hu * exp(-gamma_u_m * (1.0 - h_u * 0.8 / h_u)) # 简化为 exp(-gamma_u_m * 0.2) 亦可
  rb_u = ra_leaf_boundary(Tair, u_u)

  kh_o = 0.41 * ustar * (h - h * 0.8) / ψ # 主林冠顶部涡流扩散系数

  # h~h_u 之间的气动阻力 (介质为主林冠下部，使用 gamma_trunk)
  ra_u = h / (gamma_o_h * kh_o) * (exp(gamma_o_h * (1.0 - h_u / h)) - 1.0)

  ## 3. Ground (极近地表层)
  # K_h 从 h 衰减至 h_u (同样穿越主林冠下部，使用 gamma_trunk)
  kh_u = kh_o * exp(-gamma_o_h * (1.0 - h_u / h))

  gamma_g_h = 4.0
  # 林下层底部到地表的气动阻力 (完美的闭合公式)
  ra_g = h_u / (gamma_g_h * kh_u) * (exp(gamma_g_h) - 1.0)



  # 后处理
  ra_g = ra_g + ra_u + ra_o
  ra_g = max(120, ra_g)

  G_o_a = 1.0 / ra_o
  G_o_b = 1.0 / rb_o

  G_u_a = 1.0 / (ra_o + ra_u)
  G_u_b = 1.0 / rb_u
  return ra_o, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b
end
