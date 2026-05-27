using BEPS, Test

@testset "ParamBEPS Construction" begin
  @testset "ENF (VegType=1) with Loam (SoilType=4)" begin
    model = ParamBEPS(1, 4; N=5)

    # Test model structure
    @test model.N == 5
    @test model isa ParamBEPS{Float64}

    # Test vegetation parameters
    @test model.veg.LAI_max_o ≈ 4.5
    @test model.veg.LAI_max_u ≈ 2.4
    @test model.veg.VCmax25 ≈ 62.5
    @test model.veg.g1_w ≈ 8.0

    # Test water potential parameters (non-DBF/EBF should have higher ψ_min)
    @test model.ψ_min ≈ 33.0
    @test model.alpha ≈ 0.4

    # Test drainage and root decay
    @test model.r_drainage ≈ 0.5
    @test model.r_root_decay ≈ 0.95

    # Test soil hydraulic parameters (loam)
    @test length(model.hydraulic.θ_sat) == 5
    @test model.hydraulic.θ_sat[1] ≈ 0.463
    @test model.hydraulic.θ_res[1] ≈ 0.12
    @test model.hydraulic.b[1] ≈ 4.5

    # Test soil thermal parameters
    @test length(model.thermal.κ_dry) == 5
    @test model.thermal.κ_dry[1] ≈ 7.0
    @test model.thermal.ρ_soil[1] ≈ 1300.0
  end

  @testset "DBF (VegType=6) with Sandy Loam (SoilType=3)" begin
    model = ParamBEPS(6, 3; N=5)

    # Test DBF-specific vegetation parameters
    @test model.veg.VCmax25 ≈ 57.7
    @test model.veg.N_leaf ≈ 1.74 + 0.71

    # Test water potential parameters (DBF should have lower ψ_min)
    @test model.ψ_min ≈ 10.0
    @test model.alpha ≈ 1.5

    # Test soil hydraulic parameters (sandy loam)
    @test model.hydraulic.θ_sat[1] ≈ 0.453
    @test model.hydraulic.θ_res[1] ≈ 0.10
    @test model.hydraulic.b[1] ≈ 3.1
    @test model.hydraulic.K_sat[1] ≈ 2.592

    # Test soil thermal parameters
    @test model.thermal.κ_dry[1] ≈ 8.0
  end

  @testset "EBF (VegType=9) with Clay (SoilType=11)" begin
    model = ParamBEPS(9, 11; N=5)

    # Test EBF-specific vegetation parameters
    @test model.veg.VCmax25 ≈ 29.0

    # Test water potential parameters (EBF should have lower ψ_min like DBF)
    @test model.ψ_min ≈ 10.0
    @test model.alpha ≈ 1.5

    # Test soil hydraulic parameters (clay)
    @test model.hydraulic.θ_sat[1] ≈ 0.475
    @test model.hydraulic.θ_res[1] ≈ 0.27
    @test model.hydraulic.b[1] ≈ 7.6

    # Test soil thermal parameters
    @test model.thermal.κ_dry[1] ≈ 4.4
  end

  @testset "C4 (VegType=40) with Sand (SoilType=1)" begin
    model = ParamBEPS(40, 1; N=5)

    # Test C4-specific vegetation parameters
    @test model.veg.VCmax25 ≈ 30.0
    @test model.veg.g1_w ≈ 4.0  # C4 has different Ball-Berry coefficient

    # Test soil hydraulic parameters (sand)
    @test model.hydraulic.θ_sat[1] ≈ 0.437
    @test model.hydraulic.θ_res[1] ≈ 0.03
    @test model.hydraulic.b[1] ≈ 1.7
    @test model.hydraulic.K_sat[1] ≈ 20.88

    # Test soil thermal parameters
    @test model.thermal.κ_dry[1] ≈ 8.6
  end

  @testset "Different N layers" begin
    model3 = ParamBEPS(1, 4; N=3)
    model5 = ParamBEPS(1, 4; N=5)

    @test model3.N == 3
    @test model5.N == 5
    @test length(model3.hydraulic.θ_sat) == 3
    @test length(model5.hydraulic.θ_sat) == 5
  end

  @testset "Different Float Types" begin
    model_f64 = ParamBEPS(1, 4; N=5, FT=Float64)
    model_f32 = ParamBEPS(1, 4; N=5, FT=Float32)

    @test model_f64 isa ParamBEPS{Float64}
    @test model_f32 isa ParamBEPS{Float32}
    @test eltype(model_f64.hydraulic.θ_sat) == Float64
    @test eltype(model_f32.hydraulic.θ_sat) == Float32
  end

  @testset "Model Display" begin
    model = ParamBEPS(1, 4; N=5)

    # Test that show method works without error
    io = IOBuffer()
    show(io, model)
    output = String(take!(io))

    @test !isempty(output)
    @test occursin("ParamBEPS", output)
  end
end

@testset "InitParam_Soil Function" begin
  using BEPS: InitParam_Soil

  @testset "All Soil Types" begin
    soil_types = [
      (1, "sand"),
      (2, "loamy sand"),
      (3, "sandy loam"),
      (4, "loam"),
      (5, "silty loam"),
      (6, "sandy clay loam"),
      (7, "clay loam"),
      (8, "silty clay loam"),
      (9, "sandy clay"),
      (10, "silty clay"),
      (11, "clay")
    ]

    for (soil_type, name) in soil_types
      hydraulic, thermal = InitParam_Soil(soil_type, 5, Float64)

      @test hydraulic isa BEPS.HydraulicProfile{Float64,5}
      @test thermal isa BEPS.ThermalProfile{Float64,5}
      @test length(hydraulic.θ_sat) == 5
      @test length(thermal.κ_dry) == 5

      # All values should be positive (except ψ_sat which is negative)
      @test all(hydraulic.θ_sat .> 0)
      @test all(hydraulic.θ_res .> 0)
      @test all(hydraulic.K_sat .> 0)
      @test all(hydraulic.ψ_sat .> 0)
      @test all(hydraulic.b .> 0)
      @test all(thermal.κ_dry .> 0)
      @test all(thermal.ρ_soil .> 0)
    end
  end

  @testset "Physical Constraints" begin
    hydraulic, thermal = InitParam_Soil(4, 5, Float64)

    # θ_res < θ_sat (wilting point less than saturation)
    @test all(hydraulic.θ_res .< hydraulic.θ_sat)

    # Reasonable ranges
    @test all(0 .< hydraulic.θ_sat .< 1)
    @test all(0 .< hydraulic.θ_res .< 1)
    @test all(hydraulic.K_sat .> 0)
    @test all(hydraulic.b .> 0)
  end
end
