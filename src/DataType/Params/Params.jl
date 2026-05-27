export ParamVeg
export ParamSoilHydraulic, ParamSoilThermal, ParamSoil
export ParamBEPS, BEPSCampbell, BEPSCampbellLayers
export get_params, parameters, update!, get_opt_info
export filter_params, update_params!


using Parameters, DataFrames
import ModelParams: @metadata, @units, units, @bounds, bounds, split_bounds
import ModelParams: AbstractSoilModel, AbstractRetention, MultiLayer
import ModelParams: HydraulicProfile, ThermalProfile, ThermalBaseLayers, KvLayers
import ModelParams: get_params, parameters, update!, get_opt_info
import ModelParams: filter_params, update_params!

@bounds @with_kw mutable struct ParamVeg{FT<:AbstractFloat}
  # lc::Int = 1
  # Ω0::FT = 0.7             # // clumping_index
  has_understory::Bool = true      # 
  is_bforest::Bool = false         # broadleaf forest

  LAI_max_o::FT = 4.5 | (0.1, 7.0)      # LAI max for overstory
  LAI_max_u::FT = 2.4 | (0.1, 7.0)      # LAI max for understory

  α_canopy_vis::FT = 0.055 | (0.02, 0.15)  # canopy albedo visible
  α_canopy_nir::FT = 0.300 | (0.15, 0.50)  # canopy albedo near-infrared
  α_soil_sat::FT = 0.10 | (0.05, 0.20)     # albedo of saturated soil
  α_soil_dry::FT = 0.35 | (0.20, 0.50)     # albedo of dry soil

  # r_drainage::FT = 0.5     # ? 产流比例
  r_root_decay::FT = Cdouble(0.95) | (0.85, 0.999) # ? 根系分布衰减率, decay_rate_of_root_distribution

  Ω::FT = 0.85 | (0.3, 1.0)             # clumping_index
  z_canopy_o::FT = 1.0 | (0.1, 50.0)    # overstory canopy height [m]
  z_canopy_u::FT = 0.2 | (0.05, 5.0)    # understory canopy height [m]
  z_wind::FT = 2 | (1.0, 100.0)         # wind measurement height [m]

  g1_w::FT = 8 | (1.0, 20.0)            # Ball-Berry slope coefficient
  g0_w::FT = 0.0175 | (0.001, 0.1)      # Ball-Berry intercept for H2O

  VCmax25::FT = 89.45 | (5.0, 200.0)    # max Rubisco capacity at 25℃ [μmol m-2 s-1], global range
  # Jmax25::FT = 2.39 * 57.7 - 14.2 #

  # # coefficient reflecting the sensitivity of stomata to VPD/moderately N-stressed plants
  N_leaf::FT = 1.74 + 0.71 | (0.5, 5.0)   # leaf Nitrogen content, mean value + 1 SD [g/m2]
  slope_Vc::FT = 33.79 / 57.7 | (0.3, 1.0) # slope for Vcmax-N relationship
end


# 水力参数
@bounds @with_kw mutable struct ParamSoilHydraulic{FT<:AbstractFloat}
  θ_vfc::FT = FT(0.30) | (0.10, 0.45)   # volumetric field capacity [-]
  θ_vwp::FT = FT(0.10) | (0.02, 0.30)   # volumetric wilting point [-]

  θ_sat::FT = FT(0.45) | (0.25, 0.70)   # volumetric saturation [-]
  K_sat::FT = FT(5.0) | (0.01, 50.0)   # saturated hydraulic conductivity [cm h-1]
  ψ_sat::FT = FT(-0.5) | (-2.0, -0.01)  # matric potential at saturation [m]
  b::FT = FT(5.0) | (1.5, 15.0)    # Campbell parameter [-]
end

# 热力参数
@bounds @with_kw mutable struct ParamSoilThermal{FT<:AbstractFloat}
  κ_dry::FT = FT(0.2) | (0.05, 0.5)      # dry soil thermal conductivity [W m-1 K-1]
  ρ_soil::FT = FT(1300.0) | (800.0, 1800.0) # soil bulk density [kg m-3]
  V_SOM::FT = FT(0.02) | (0.0, 0.3)      # organic matter volume fraction [-]
end

@with_kw mutable struct ParamSoil{FT<:AbstractFloat}
  hydraulic::ParamSoilHydraulic{FT} = ParamSoilHydraulic{FT}()
  thermal::ParamSoilThermal{FT} = ParamSoilThermal{FT}()
end


include("GlobalData.jl")
include("Param_Init.jl")

include("../BEPS_Param.jl")

include("deprecated/VegHelper.jl")
include("deprecated/ReadParamVeg.jl")
include("deprecated/Init_Soil_Parameters.jl")

include("ParamPhoto.jl")
