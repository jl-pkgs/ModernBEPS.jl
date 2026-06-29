@testset "aerodynamic_conductance" begin
  canopyh_o = 2.0
  canopyh_u = 0.2
  height_wind_sp = 10.0   # 塔测，远在冠层之上，不触发 safe_wind_ref 调整
  clumping = 0.8
  Ta = 20.0
  wind_sp = 2.0
  GH_o = 100.0
  pai_o = 4.0
  pai_u = 2.0

  ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
    aerodynamic_conductance_jl(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
      pai_o, pai_u)
  r1 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

  ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u =
    clang.aerodynamic_conductance_c(canopyh_o, canopyh_u, height_wind_sp, clumping, Ta, wind_sp, GH_o,
      pai_o, pai_u)
  r2 = (; ra_o, ra_u, ra_g, Ga_o, Gb_o, Ga_u, Gb_u)

  @test maximum(abs.(values(r1) .- values(r2))) <= 1e-8
end
