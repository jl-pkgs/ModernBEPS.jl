# const RTIMES = 24.0 # 呼吸作用系数转换常数 (day -> hour)

"""
    InitParam_Veg(lc::AbstractString; FT=Float64)

读取 JSON 配置文件并返回 ParamVeg 结构体。
"""
function InitParam_Veg(VegType::AbstractString="ENF"; FT=Float64)
  veg_data = JSON.parsefile(PATH_VEG)
  gen_data = JSON.parsefile(PATH_GEN)

  VegCode = find_VegType(VegType)
  v = veg_data[veg_param_key(VegType)]

  return ParamVeg{FT}(
    has_understory = !(VegCode == 25 || VegCode == 40),
    is_bforest    = VegCode == 6 || VegCode == 9,
    LAI_max_o    = FT(v["LAI_max_o"]),
    LAI_max_u    = FT(v["LAI_max_u"]),
    α_canopy_vis = FT(v["albedo_canopy_vis"]),
    α_canopy_nir = FT(v["albedo_canopy_nir"]),
    α_soil_sat   = FT(gen_data["albedo_saturated_soil"]),
    α_soil_dry   = FT(gen_data["albedo_dry_soil"]),
    z_canopy_o   = FT(v["z_canopy_o"]),
    z_canopy_u   = FT(v["z_canopy_u"]),
    z_wind       = FT(gen_data["the_height_to_measure_wind_speed"]),
    g1_w         = FT(v["g1_w"]),
    g0_w         = FT(gen_data["intercept_for_H2O_ball_berry"]),
    VCmax25      = FT(v["VCmax25"]),
    N_leaf       = FT(v["N_leaf"]),
    slope_Vc     = FT(v["slope_Vc"])
  )
end

InitParam_Veg(VegType::Integer; FT=Float64) = InitParam_Veg(veg_name(VegType); FT)


"""
    InitParam_Soil(SoilType::AbstractString, N::Int, FT::Type)

Initialize soil hydraulic and thermal parameters.
SoilType: 1=sand, 2=loamy sand, 3=sandy loam, 4=loam, 5=silty loam,
          6=sandy clay loam, 7=clay loam, 8=silty clay loam,
          9=sandy clay, 10=silty clay, 11=clay
"""
function InitParam_Soil(SoilType::AbstractString, N::Int, FT::Type)
  idx = find_SoilType(SoilType)
  InitParam_Soil(idx, N, FT)
end

function InitParam_Soil(SoilType::Integer, N::Int, FT::Type)
  idx = find_SoilType(SoilType)
  p = SOIL_PARAMS[idx]

  b = _fit_layers(FT, p.b, N)            # [-], campbell's b parameter

  K_sat = _fit_layers(FT, p.K_sat, N)    # [cm h-1]
  θ_sat = fill(FT(p.θ_sat), N) # [%]
  # θ_vfc = fill(FT(p.θ_vfc), n) # [%]
  θ_vwp = fill(FT(p.θ_vwp), N) # [%]
  ψ_sat = _fit_layers(FT, p.ψ_sat, N)    # [m], positive suction at saturation (Campbell 1974 convention)

  SOIL_THERMAL_DENSITY = [1300.0, 1500.0, 1517.0, 1517.0, 1517.0] # [kg m-3]
  SOIL_ORGANIC_MATTER = [0.05, 0.02, 0.01, 0.01, 0.003]           # volume fraction, 0-1

  κ_dry = fill(FT(p.κ_dry), N) # [W m-1 K-1]
  ρ_soil = _fit_layers(FT, SOIL_THERMAL_DENSITY, N) # [kg m-3]
  V_SOM = _fit_layers(FT, SOIL_ORGANIC_MATTER, N)   # [volume fraction], 0-1

  dz = _default_dz(FT, N)
  profile = BEPSCampbellLayers{FT,N}(; θ_vwp, θ_sat, Ksat=K_sat, ψ_sat, b)
  kv = KvLayers{FT,N}(; kv=K_sat)
  hydraulic = HydraulicProfile{FT,N}(profile, kv, FT.(100 .* dz))
  thermal = ThermalProfile{FT,N}(ThermalBaseLayers{FT,N}(; κ_dry, ρ_soil, V_SOM))
  return hydraulic, thermal
end

# if VegType == 6 || VegType == 9 # DBF or EBF, low constaint threshold
#   p.ψ_min = 10.0 # ψ_min
#   p.alpha = 1.5
# else
#   p.ψ_min = 33.0 # ψ_min
#   p.alpha = 0.4
# end
