"""
    aerodynamic_conductance_jl(canopy_height_o::FT, canopy_height_u::FT,
        z_wind::FT, clumping::FT,
        Tair::FT, wind_sp::FT, SH_o_p::FT, lai_o::FT, lai_u::FT=0.0)

# Examples
```julia
canopyh_o = 2.0
canopyh_u = 0.2
height_wind_sp = 2.0
clumping = 0.8
Ta = 20.0
wind_sp = 2.
GH_o = 100.0
lai_o = 4.0
lai_u = 2.0

ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
  aerodynamic_conductance_jl(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
    lai_o, lai_u)
r1 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)
```
"""
function aerodynamic_conductance_jl(canopy_height_o::FT, canopy_height_u::FT,
  z_wind::FT, clumping::FT,
  Tair::FT, wind_sp::FT, SH_o_p::FT, lai_o::FT, lai_u::FT=0.0)

  # beta::FT = 0.5 # Bowen's ratio
  k::FT = 0.4   # von Karman's constant
  cp::FT = 1010 # specific heat of air (J/kg/K)
  density_air::FT = 1.225 # density of air at 15 C (kg/m3)
  gg::FT = 9.8 # gravitational acceleration (m/s2)
  n::FT = 5.0
  nu_lower::FT = (13.3 + Tair * 0.07) / 1000000  # viscosity (cm2/s)
  alfaw::FT = (18.9 + Tair * 0.07) / 1000000

  wind_sp = max(wind_sp, 0.02)
  @fastmath if !(isfinite(wind_sp) && wind_sp > 0)
    G_o_a = 1 / 200.0
    G_o_b = 1 / 200.0
    G_u_a = 1 / 400.0
    G_u_b = 1 / 200.0
    ra_g = 600.0
    ra_u = 200.0
    ra_o = 200.0
  else
    d::FT = 0.8 * canopy_height_o  # displacement height (m)
    z0::FT = 0.08 * canopy_height_o # roughness length (m)
    # 确保参考高度在 d+z0 之上（再分析10m风 + 高冠层保护）
    z_eff, _ = safe_wind_ref(z_wind, canopy_height_o; d_frac=FT(0.8), z0m_frac=FT(0.08))

    ustar::FT = wind_sp * k / log((z_eff - d) / z0) # friction velocity (m/s)
    L::FT = -(k * gg * SH_o_p) / (density_air * cp * (Tair + 273.3) * ustar^3)
    L = max(-2.0, L)

    ra_o::FT = 1 / (k * ustar) * (log((z_eff - d) / z0) + (n * (z_eff - d) * L))
    ra_o = clamp(ra_o, 2, 100)
    
    if L > 0
      ψ = 1 + 5 * (z_eff - d) * L
    else
      ψ = (1 - 16 * (z_eff - d) * L)^(-0.5)
    end
    ψ = min(10.0, ψ)

    #******** Leaf boundary layer resistance ******************/
    # Wind speed at tree top */
    uh::FT = 1.1 * ustar / k  # wind speed at height h
    Le = lai_o * clumping
    gamma = (0.167 + 0.179 * uh) * Le^(1.0 / 3.0)

    # Wind speed at d, taking as the mean wind speed inside a stand */
    ud = uh * exp(-gamma * (1 - d / canopy_height_o))

    Nu = cal_Nu(ud, nu_lower)
    rb_o = min(40, 0.5 * 0.1 / (alfaw * Nu)) # leaf boundary resistance

    # uf = ustar
    G_o_a = 1.0 / ra_o
    G_o_b = 1.0 / rb_o

    gamma = 0.1 + lai_o^0.75
    kh_o = 0.41 * ustar * (canopy_height_o - canopy_height_o * 0.8) / ψ

    # wind speed at the zero displacement of canopy
    un_d = uh * exp(-gamma * (1 - canopy_height_u * 0.8 / canopy_height_o))
    # # wind speed at the zero displacement of canopy
    # un_t = uh * exp(-gamma * (1 - canopy_height_u / canopy_height_o))
    Nu = cal_Nu(un_d, nu_lower)
    rb_u = min(40, 0.5 * 0.1 / (alfaw * Nu)) # leaf boundary resistance
    G_u_b = 1.0 / rb_u

    ra_u = canopy_height_o / (gamma * kh_o) * (exp(gamma * (1 - canopy_height_u / canopy_height_o)) - 1)
    G_u_a = 1.0 / (ra_o + ra_u)

    gamma = 4.0
    # kh_u = kh_o * exp(-gamma * (1 - canopy_height_u / canopy_height_o))
    ra_g = canopy_height_o / (gamma * kh_o) * (exp(gamma * (1 - 0 / canopy_height_o)) -
                                               exp(gamma * (1 - canopy_height_u / canopy_height_o)))
    ra_g = ra_g + ra_u + ra_o
    ra_g = max(120, ra_g)
  end

  return ra_o, ra_u, ra_g, G_o_a, G_o_b, G_u_a, G_u_b
end

# function windProfile_factor(canopy_height_u, canopy_height_o, gamma, k=1.0)
#   exp(gamma * (1 - canopy_height_u * k / canopy_height_o))
# end

function cal_Nu(u::FT, nu_lower::FT)::FT
  # lw::T = 0.3 # leaf characteristic width =0.3 for BS
  # sigma::T = 5 # shelter factor =5 for BS
  # Re = (ud * lw / sigma) / nu_lower
  Re::FT = (u * 0.1) / nu_lower # Reynold's number
  Nu::FT = 1.0 * Re^0.5 # Nusselt number
  Nu
end
