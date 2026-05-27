using BEPS, Test


@testset "Model Parameters" begin
  # xs = ParamSoilHydraulicLayers{Float64,4}()
  model = ParamBEPS{Float64}(; N=5)
  params = parameters(model)
  display(model)

  update!(model, [[:r_drainage]], [0.4]; params)
  @test model.r_drainage == 0.4

  p_hydraulic = filter_params(model, :hydraulic)
  update_params!(model, [[:hydraulic, :profile, :b, 4]], [4.0]; params=p_hydraulic)
  @test model.:hydraulic.b[4] == 4.0

  @test all(path -> path[1] === :hydraulic, p_hydraulic.path)
end

# @testset "ParamSoilHydraulicLayers" begin
#   x = ParamSoilHydraulic{Float64}()
#   xs = ParamSoilHydraulicLayers{Float64,4}()
#   length(get_bounds(x)) == 6
#   length(get_bounds(xs)) == 24
# end
