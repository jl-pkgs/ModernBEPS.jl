using BEPS, Test


@testset "Model Parameters" begin
  model = ParamBEPS{Float64}(; N=5)
  params = parameters(model)
  display(model)

  paths = [
    [:r_drainage],
    [:hydraulic, :profile, :b, 4]
  ]
  values = [0.4, 4.0]

  update!(model, paths, values; params)
  @test model.r_drainage == 0.4
  @test model.:hydraulic.b[4] == 4.0

  p_hydraulic = filter_params(model, :hydraulic)
  @test all(path -> path[1] === :hydraulic, p_hydraulic.path)
  @test size(p_hydraulic, 1) == 25
end
