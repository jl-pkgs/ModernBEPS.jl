using BEPS, Test
import BEPS: Params2Soil!, Soil2Params!, State2Soil!, UpdateRootFraction!, Init_Soil_T_θ!
import BEPS: UpdateThermal_κ, UpdateThermal_Cv, Update_ice_ratio, UpdateHeatFlux, UpdateSoilMoisture, Root_Water_Uptake, Init_Soil_var

kstep = 360.0

function Base.getindex(x::Union{Soil,Soil_c}, i::Integer)
  names = fieldnames(typeof(x))
  getfield(x, names[i])
end

function is_soil_equal(p_jl, p_c; tol=1e-7, verbose=false)
  names_c = fieldnames(typeof(p_c))
  names_jl = fieldnames(typeof(p_jl))
  names_skip = [:θb, :ψb, :θ_vfc]

  for i in eachindex(names_c)
    name_c = names_c[i]
    name_jl = names_jl[i]
    name_jl in names_skip && continue

    x_jl = getfield(p_jl, i)
    x_c = getfield(p_c, i)

    if name_jl == :K_sat
      # C version uses m/s, Julia uses cm/h
      x_jl = x_jl ./ 360000.0
    elseif name_jl in (:Kmid, :Kavg)
      x_jl = x_jl ./ 360000.0
    elseif name_jl == :ψ
      x_jl = .-x_jl ./ 100.0
    elseif name_jl == :ψ_sat
      x_jl = .-x_jl ./ 100.0
    elseif name_jl == :ψ_min
      x_jl = x_jl ./ 100.0
    end

    if verbose
      ## 变量名可能不同
      color = name_c == name_jl ? :black : :red
      printstyled("C and Julia: $name_c, $name_jl\n"; color)
    end
    @test maximum(abs.(x_c .- x_jl)) <= tol
  end
end


# Init_Soil_Parameters(par.landcover, par.soil_type, parameter[28], p_soil)
@testset "Init_Soil_Parameters" begin
  p_jl = Soil()
  p_c = Soil_c()
  # stxt = 7: fix const κ_dry
  # stxt = 6: fix K_sat
  for stxt = [1, 2, 3, 4, 5, 8, 9, 10, 11] # skip 6, 7
    Init_Soil_Parameters(p_jl, 1, stxt, 0.1)
    clang.Init_Soil_Parameters(p_c, 1, stxt, 0.1)
    is_soil_equal(p_jl, p_c; tol=1e-8, verbose=false)
  end
end

function init_soil()
  kw = (lon=120.5, lat=30.5,
    VegType=25, SoilType=8,
    clumping=0.85,
    Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
  )
  param = ReadParamVeg(kw.VegType)      # n = 48

  rad = 100.0
  Tair = -4.0
  hum = 7.0
  pre = 0.0
  wind = 2.0
  meteo = Met(rad, 0.0, Tair, hum, pre, wind)

  Tsoil_p = collect(4.0:-1:-5.0)
  Tsoil_c = collect(-5.0:1:4.0)
  ice_ratio = fill(0.6, 10)

  p_jl = Soil(; Tsoil_p, Tsoil_c, ice_ratio)
  p_c = Soil_c(;
    Tsoil_p=Tuple(Tsoil_p),
    Tsoil_c=Tuple(Tsoil_c),
    ice_ratio=Tuple(ice_ratio))

  # var_jl = zeros(41)
  state_jl = StateBEPS()
  state_c = zeros(41)
  Ta = meteo.Tair

  r_drainage = param[27]
  r_root_decay = param[28]

  Init_Soil_var(p_jl, state_jl, Ta; r_drainage, r_root_decay, kw...)
  Init_Soil_var(p_c, state_c, Ta; r_drainage, r_root_decay, kw...)
  p_jl, p_c
end


@testset "UpdateHeatFlux" begin
  p_jl, p_c = init_soil()

  funs = [
    # UpdateThermal_Cv,
    # Update_ice_ratio,
    UpdateThermal_κ,
    soil_water_factor_v2]

  for fun in funs
    fun(p_jl)
    fun(p_c)
  end
  
  UpdateHeatFlux(p_jl, 20.0, 3600.0)
  UpdateHeatFlux(p_c, 20.0, 3600.0)
  is_soil_equal(p_jl, p_c; tol=1e-7, verbose=true)
end


@testset "Init_Soil_var" begin
  p_jl, p_c = init_soil()
  funs = [
    UpdateThermal_Cv,
    Update_ice_ratio,
    UpdateThermal_κ,
    soil_water_factor_v2]

  for fun in funs
    fun(p_jl)
    fun(p_c)
  end

  UpdateSoilMoisture(p_jl, kstep)
  UpdateSoilMoisture(p_c, kstep)

  UpdateHeatFlux(p_jl, 20.0, 3600.0)
  UpdateHeatFlux(p_c, 20.0, 3600.0)
  is_soil_equal(p_jl, p_c; tol=1e-7, verbose=true)
end


# 测试新旧 API 兼容性：SoilState + ParamBEPS vs Soil
@testset "StateBEPS + ParamBEPS API 兼容性" begin
  # 创建 ParamBEPS 参数
  ps = ParamBEPS(25, 8)  # VegType=25, SoilType=8

  # 创建旧版 Soil 结构
  soil = Soil()
  Params2Soil!(soil, ps)
  Init_Soil_T_θ!(soil, 2.2, 10.0, 0.4, 0.0)
  UpdateRootFraction!(soil)

  # 创建新版 StateBEPS
  st = StateBEPS()
  st.n_layer = Cint(ps.N)
  st.dz[1:ps.N] .= ps.dz
  Init_Soil_T_θ!(st, 2.2, 10.0, 0.4, 0.0)
  UpdateRootFraction!(st, ps)

  # 验证两种 API 结果一致
  @test st.z_snow ≈ soil.z_snow
  @test st.f_root[1:5] ≈ soil.f_root[1:5]
  @test st.Tsoil_c[1:5] ≈ soil.Tsoil_c[1:5]
  @test st.Tsoil_p[1:5] ≈ soil.Tsoil_p[1:5]
  @test st.θ[1:5] ≈ soil.θ[1:5]
  @test st.ice_ratio[1:5] ≈ soil.ice_ratio[1:5]

  # 测试 soil_water_factor_v2
  soil_water_factor_v2(soil)
  soil_water_factor_v2(st, ps)

  @test st.f_soilwater ≈ soil.f_soilwater
  @test st.f_stress[1:5] ≈ soil.f_stress[1:5]
  @test st.w_norm[1:5] ≈ soil.w_norm[1:5]

  # 测试 UpdateThermal_κ
  UpdateThermal_κ(soil)
  UpdateThermal_κ(st, ps)
  @test st.κ[1:5] ≈ soil.κ[1:5]

  # 测试 UpdateThermal_Cv
  UpdateThermal_Cv(soil)
  UpdateThermal_Cv(st, ps)
  @test st.Cv[1:5] ≈ soil.Cv[1:5]

  # 测试 UpdateHeatFlux (implicit Update_ice_ratio)
  UpdateHeatFlux(soil, 20.0, 3600.0)
  UpdateHeatFlux(st, 20.0, 3600.0)
  @test st.Tsoil_c[1:5] ≈ soil.Tsoil_c[1:5]
  @test st.Tsoil_p[1:5] ≈ soil.Tsoil_p[1:5]
  @test st.G[1:5] ≈ soil.G[1:5]
  @test st.ice_ratio[1:5] ≈ soil.ice_ratio[1:5]

  # 测试 UpdateSoilMoisture
  UpdateSoilMoisture(soil, 3600.0)
  UpdateSoilMoisture(st, ps, 3600.0)
  @test st.θ[1:5] ≈ soil.θ[1:5]
  @test st.z_water ≈ soil.z_water
  @test st.r_waterflow[1:5] ≈ soil.r_waterflow[1:5]
  @test st.ψ[1:5] ≈ soil.ψ[1:5]
  @test all(st.ψ[1:5] .<= 0.0)
  @test all(st.Kmid[1:5] .>= 0.0)
  @test all(st.Kavg[1:4] .>= 0.0)

  # 测试 Root_Water_Uptake
  Root_Water_Uptake(soil, 1.0, 2.0, 0.5)
  Root_Water_Uptake(st, 1.0, 2.0, 0.5)
  @test st.ETi[1:5] ≈ soil.ETi[1:5]
end


@testset "Soil → StateBEPS 转换" begin
  # 创建并初始化 Soil
  ps = ParamBEPS(25, 8)
  soil = Soil()
  Params2Soil!(soil, ps)
  Init_Soil_T_θ!(soil, 2.2, 10.0, 0.4, 0.0)
  UpdateRootFraction!(soil)
  soil_water_factor_v2(soil)

  # 从 Soil 构造 StateBEPS
  st = StateBEPS(soil)

  # 验证所有状态变量一致
  @test st.n_layer == soil.n_layer
  @test st.z_water ≈ soil.z_water
  @test st.z_snow ≈ soil.z_snow
  @test st.f_soilwater ≈ soil.f_soilwater
  @test st.θ[1:5] ≈ soil.θ[1:5]
  @test st.Tsoil_c[1:5] ≈ soil.Tsoil_c[1:5]
  @test st.ice_ratio[1:5] ≈ soil.ice_ratio[1:5]
  @test st.f_root[1:5] ≈ soil.f_root[1:5]

  # 测试 State2Soil! 回写
  st.θ[1] = 0.5
  st.Tsoil_c[1] = 15.0
  State2Soil!(soil, st)
  @test soil.θ[1] ≈ 0.5
  @test soil.Tsoil_c[1] ≈ 15.0
end


@testset "Soil2Params! 逆函数测试" begin
  # Initialize with default
  ps = ParamBEPS(25, 8)
  soil = Soil()
  Params2Soil!(soil, ps)
  
  # Modify soil parameters
  soil.r_drainage = 0.99
  soil.ψ_min = 9900.0
  soil.alpha = 9.9
  # soil.θ_vfc[1] = 0.88
  soil.K_sat[1] = 0.77
  
  # Sync back to params
  Soil2Params!(ps, soil)
  
  @test ps.r_drainage ≈ 0.99
  @test ps.ψ_min ≈ 9900.0
  @test ps.alpha ≈ 9.9
  # @test ps.hydraulic.θ_vfc[1] ≈ 0.88
  @test ps.hydraulic.K_sat[1] ≈ 0.77
end
