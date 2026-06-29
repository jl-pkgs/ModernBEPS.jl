using Test, BEPS

## V1: 独立模块，保留旧版用于对比（已被 BEPS 模块替换为 V2）
include("diagnose_ra.jl")

module AC_V1
using Parameters: @with_kw
using UnPack: @pack!
const FT = Float64
include(normpath(joinpath(@__DIR__, "../../src/DataType/AeroConsts.jl")))
include(normpath(joinpath(@__DIR__, "../../src/aerodynamic_conductance.jl")))
end

## V2: 独立模块，不依赖 BEPS 导出，方便物理验证
module AC_V2
using Parameters: @with_kw
using UnPack: @pack!
const FT = Float64
include(normpath(joinpath(@__DIR__, "../../src/DataType/AeroConsts.jl")))
include(normpath(joinpath(@__DIR__, "../../src/aerodynamic_conductance_V2.jl")))
end

# ── V1 vs C（仅 Windows 有效） ──────────────────────────────────────

@testset "aerodynamic_conductance V1 vs C" begin
  canopyh_o = 2.0
  canopyh_u = 0.2
  height_wind_sp = 10.0   # 远在冠层之上，不触发 safe_wind_ref，确保与 C 库可比
  clumping = 0.8
  Ta = 20.0
  wind_sp = 2.0
  GH_o = 100.0
  pai_o = 4.0
  pai_u = 2.0

  ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
    AC_V1.aerodynamic_conductance_jl(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
      pai_o, pai_u)
  r1 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

  ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
    clang.aerodynamic_conductance_c(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
      pai_o, pai_u)
  r2 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

  @test maximum(abs.(values(r1) .- values(r2))) <= 1e-8
end

# ── V2 物理正确性测试（跨平台，不需要 C 库）──────────────────────────
@testset "aerodynamic_conductance V2 物理正确性" begin
  # 典型高大森林参数
  h = 20.0;
  h_u = 2.0;
  z_wind = 40.0
  Ω = 0.8;
  Ta = 20.0;
  u = 3.0
  lai_o = 4.0;
  lai_u = 2.0

  # 1. 零风速 → 返回默认值
  ra_o0, ra_u0, ra_g0, Ga_o0, Gb_o0, Ga_u0, Gb_u0 =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, 0.0, 100.0, lai_o, lai_u)
  @test ra_o0 == 0.0
  @test ra_g0 == 300.0
  @test Ga_o0 ≈ 1 / 200.0
  @test Ga_u0 ≈ 1 / 200.0

  # 2. 近中性（H≈0）→ ra_o 接近解析解
  #    u* ≈ u·k/ln((z-d)/z0m), ra_o_neutral ≈ 1/(k·u*)·ln((z-d)/z0h)
  k = 0.4;
  d = 0.8h;
  z0m = 0.08h;
  z0h = 0.1z0m
  z_eff, _ = AC_V2.safe_wind_ref(z_wind, h)
  ustar_neutral = u * k / log((z_eff - d) / z0m)
  ra_o_expected = 1 / (k * ustar_neutral) * log((z_eff - d) / z0h)

  ra_o_n, ra_u_n, ra_g_n, Ga_o_n, Gb_o_n, Ga_u_n, Gb_u_n =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, 1e-6, lai_o, lai_u)
  @test ra_o_n ≈ ra_o_expected rtol = 0.01   # 1% 偏差

  # 3. 物理合理性：稳定 (H<0) > 中性 > 不稳定 (H>0) → ra_o 有序
  ra_o_stable = first(AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, -200.0, lai_o, lai_u))
  ra_o_unstable = first(AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, 200.0, lai_o, lai_u))
  @test ra_o_stable >= ra_o_n          # 稳定层结阻力更大
  @test ra_o_unstable <= ra_o_n        # 不稳定层结阻力更小

  # 4. 冠层各层关系：Ga_o > Ga_u（上层冠层更靠近参考高度，阻力更小）
  ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, 100.0, lai_o, lai_u)
  @test Ga_o > Ga_u                    # 上层传输导度更大
  @test ra_g >= 120                    # 强制最小值
  @test ra_u > 0
  @test Gb_o > 0
  @test Gb_u > 0

  # 5. 叶片边界层导度有限（非零非无穷）
  @test 0 < Gb_o < 1.0
  @test 0 < Gb_u < 1.0
end

# ── V1 vs V2 差异对比（记录已知差异） ───────────────────────────────
@testset "aerodynamic_conductance V1 vs V2 差异对比" begin
  # 使用典型森林参数
  h = 20.0;
  h_u = 2.0;
  z_wind = 40.0;
  Ω = 0.8
  Ta = 20.0;
  u = 3.0;
  H = 100.0;
  lai_o = 4.0;
  lai_u = 2.0

  ra_o1, ra_u1, ra_g1, Ga_o1, Gb_o1, Ga_u1, Gb_u1 =
    AC_V1.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, H, lai_o, lai_u)
  ra_o2, ra_u2, ra_g2, Ga_o2, Gb_o2, Ga_u2, Gb_u2 =
    AC_V2.aerodynamic_conductance_jl(h, h_u, z_wind, Ω, Ta, u, H, lai_o, lai_u)

  # V2 使用正确 MOST 公式，结果与 V1 存在合理差异
  # 两者均应在合理物理范围内
  for ra in (ra_o1, ra_o2)
    @test 2.0 <= ra <= 500.0
  end
  for ra in (ra_u1, ra_u2)
    @test ra > 0
  end
  @test ra_g1 >= 120
  @test ra_g2 >= 120

  # V2 改进点：使用 z0h 导致 ra_o 偏大（z0h < z0m，对数项更大）
  @test ra_o2 > ra_o1
end

# ── diagnose_ra 单元测试 ────────────────────────────────────────────
@testset "ra_from_flux" begin
  ρₐ = 1.292;
  cp = 1010.0

  # 典型白天条件
  r_H = ra_from_flux(200.0, 25.0, 20.0)
  @test r_H ≈ ρₐ * cp * 5.0 / 200.0 rtol = 1e-10

  # 冷冠层（SH < 0, Tc < Ta）→ r_H > 0（阻力始终为正，两者同号）
  r_H_neg = ra_from_flux(-50.0, 15.0, 20.0)
  @test r_H_neg > 0

  # Tc == Ta → r_H = 0（无温度梯度）
  @test ra_from_flux(50.0, 20.0, 20.0) ≈ 0.0

  # |SH| < 1 → Inf（数值病态，过滤）
  @test ra_from_flux(0.5, 25.0, 20.0) == Inf
  @test ra_from_flux(-0.5, 20.0, 25.0) == Inf

  # AbstractFloat 约束：Float32 输入应正常工作
  r_H_f32 = ra_from_flux(200f0, 25f0, 20f0)
  @test r_H_f32 isa Float32
end

@testset "Gc_penman" begin
  # 典型夏日午后（草地）
  Ta = 25.0;
  RH = 60.0;
  Rn = 400.0;
  G_soil = 40.0;
  LE = 200.0;
  Ga = 0.05

  Gc = Gc_penman(LE, Ga, Rn, G_soil, Ta, RH)
  @test Gc > 0
  @test Gc < 0.1   # 典型 Gc 范围 [0.001, 0.05] m/s

  # 分母趋零 → NaN（能量平衡极端条件）
  @test isnan(Gc_penman(0.0, 0.0, 0.0, 0.0, Ta, RH))
end

@testset "stability_class" begin
  @test stability_class(-1.0) == :unstable
  @test stability_class(-0.5001) == :unstable  # 刚越过下边界
  @test stability_class(-0.5) == :neutral   # 边界 -0.5：< -0.5 才是不稳定，-0.5 本身为中性
  @test stability_class(-0.6) == :unstable
  @test stability_class(-0.4) == :neutral
  @test stability_class(0.0) == :neutral
  @test stability_class(0.4) == :neutral
  @test stability_class(0.5) == :neutral   # 边界 0.5：> 0.5 才是稳定，0.5 本身为中性
  @test stability_class(0.5001) == :stable    # 刚越过上边界
  @test stability_class(0.6) == :stable
  @test stability_class(2.0) == :stable
end
