@with_kw mutable struct Flux <: AbstractFlux
  # Component gpp_* fields are instantaneous rates on ground area.
  # Convert to hourly [gC m-2] with `x * 12 * step * 1e-6`.
  gpp_o_sunlit::Cdouble = 0.0  # [μmol m-2 s-1]
  gpp_u_sunlit::Cdouble = 0.0  # [μmol m-2 s-1]
  gpp_o_shaded::Cdouble = 0.0  # [μmol m-2 s-1]
  gpp_u_shaded::Cdouble = 0.0  # [μmol m-2 s-1], 若需[gC m-2]，需[12 * step * 1e-6]

  plant_resp::Cdouble = 0.0    # [gC m-2]
  npp_o::Cdouble = 0.0         # [gC m-2]
  npp_u::Cdouble = 0.0         # [gC m-2]
  GPP::Cdouble = 0.0           # [gC m-2], hourly total; = Σgpp * 12 * step * 1e-6
  NPP::Cdouble = 0.0           # [gC m-2]
  NEP::Cdouble = 0.0           # [gC m-2]
  soil_resp::Cdouble = 0.0     # [gC m-2]
  Net_Rad::Cdouble = 0.0       # [W m-2]

  SH::Cdouble = 0.0            # [W m-2], sensible heat flux
  LH::Cdouble = 0.0            # [W m-2], latent heat flux
  Trans::Cdouble = 0.0         # [mm], hourly total
  Evap::Cdouble = 0.0          # [mm], hourly total

  z_water::Cdouble = 0.0       # [m], ponded water depth
  z_snow::Cdouble = 0.0        # [m], snow depth
  ρ_snow::Cdouble = 0.0        # [kg m-3], snow density
  inf::Cdouble = 0.0           # [cm h-1], infiltration rate
end


@with_kw mutable struct ETFlux <: AbstractFlux
  Trans_o::Cdouble = 0.0       # [kg m-2 s-1], overstory transpiration rate
  Trans_u::Cdouble = 0.0       # [kg m-2 s-1], understory transpiration rate
  Eil_o::Cdouble = 0.0         # [kg m-2 s-1], Ei of liquid (overstory)
  Eil_u::Cdouble = 0.0         # [kg m-2 s-1], Ei of liquid (understory)
  EiS_o::Cdouble = 0.0         # [kg m-2 s-1], Ei of solid (overstory)
  EiS_u::Cdouble = 0.0         # [kg m-2 s-1], Ei of solid (understory)
  Evap_soil::Cdouble = 0.0     # [kg m-2 s-1], soil evaporation rate
  Evap_SW::Cdouble = 0.0       # [kg m-2 s-1], evaporation from water pond
  Evap_SS::Cdouble = 0.0       # [kg m-2 s-1], evaporation from snow pack
  Qhc_o::Cdouble = 0.0         # [W m-2], sensible heat from overstory canopy
  Qhc_u::Cdouble = 0.0         # [W m-2], sensible heat from understory canopy
  Qhg::Cdouble = 0.0           # [W m-2], sensible heat from ground

  # Result part (hourly totals or instantaneous)
  Trans::Cdouble = 0.0         # [mm], hourly total transpiration
  Evap::Cdouble = 0.0          # [mm], hourly total evaporation
  SH::Cdouble = 0.0            # [W m-2], total sensible heat
  LH::Cdouble = 0.0            # [W m-2], total latent heat
end


function update_ET!(x::ETFlux, mid_res::Flux, Ta)
  Lv_liquid = (2.501 - 0.00237 * Ta) * 1000000  # The latent heat of water vaporization in j/kg

  x.Trans = (x.Trans_o + x.Trans_u) * step

  x.Evap = (x.Eil_o + x.Eil_u +
            x.EiS_o + x.EiS_u +
            x.Evap_soil +
            x.Evap_SW +
            x.Evap_SS) * step
  x.LH = Lv_liquid * (x.Trans_o + x.Trans_u + x.Eil_o + +
                      x.Eil_u + x.Evap_soil + x.Evap_SW) +
         Lv_solid * (x.EiS_o + x.EiS_u + x.Evap_SS)

  x.SH = (x.Qhc_o + x.Qhc_u + x.Qhg)

  # fill values to res
  mid_res.Trans = x.Trans
  mid_res.Evap = x.Evap
  mid_res.LH = x.LH
  mid_res.SH = x.SH
end

@with_kw mutable struct Radiation <: AbstractFlux
  Rs_o_df::Cdouble = 0.0
  Rs_u_df::Cdouble = 0.0

  Rs_o_dir::Cdouble = 0.0
  Rs_u_dir::Cdouble = 0.0

  Rns_o_df::Cdouble = 0.0
  Rns_u_df::Cdouble = 0.0
  Rns_g_df::Cdouble = 0.0

  Rns_o_dir::Cdouble = 0.0
  Rns_u_dir::Cdouble = 0.0
  Rns_g_dir::Cdouble = 0.0

  Rs_df::Cdouble = 0.0
  Rs_dir::Cdouble = 0.0
end

@with_kw mutable struct Cpools
  Ccd::NTuple{3,Cdouble} = nzero(3)
  Cssd::NTuple{3,Cdouble} = nzero(3)
  Csmd::NTuple{3,Cdouble} = nzero(3)
  Cfsd::NTuple{3,Cdouble} = nzero(3)
  Cfmd::NTuple{3,Cdouble} = nzero(3)
  Csm::NTuple{3,Cdouble} = nzero(3)
  Cm::NTuple{3,Cdouble} = nzero(3)
  Cs::NTuple{3,Cdouble} = nzero(3)
  Cp::NTuple{3,Cdouble} = nzero(3)
end


# ── 时间序列容器（用 @DefFluxSeries 自动展开为 Vector{FT} 字段）──────────────
@DefFluxSeries FluxSeries = Flux
@DefFluxSeries ETSeries = ETFlux
@DefFluxSeries RadSeries = Radiation

export FluxSeries, ETSeries, RadSeries
