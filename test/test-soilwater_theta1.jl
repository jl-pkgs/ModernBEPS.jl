using BEPS, Test, Dates, DataFrames

@testset "θ1-driven soilwater calibration" begin
  forcing_all = deserialize(path_proj("data/p1_forcing"))
  ntime = 24
  forcing = MetSeries(;
    ntime,
    Rs=forcing_all.Rs[1:ntime],
    Rln_in=forcing_all.Rln_in[1:ntime],
    Tair=forcing_all.Tair[1:ntime],
    RH=forcing_all.RH[1:ntime],
    Prcp=forcing_all.Prcp[1:ntime],
    Uz=forcing_all.Uz[1:ntime])
  dates = collect(DateTime(2010):Hour(1):DateTime(2010, 1, 1, 23))

  model = ParamBEPS(25, 8)
  state = InitState0(model, forcing)

  θ1_obs = fill(model.hydraulic.θ_sat[1] * 1.2, ntime)
  df = simulate_soilwater(forcing, dates; ps=model, state, θ1_obs)
  @test size(df, 1) == ntime
  @test all(df[!, :θ1] .≈ model.hydraulic.θ_sat[1])

  depths_SM = [0.15, 0.30, 0.60, 1.0]
  vars_SM = map(d -> Symbol("SM_$(Int(d * 100))cm"), depths_SM)
  SM_sim_mat = hcat([df[!, Symbol("θ$j")] for j in 1:5]...)
  SM_obs = DataFrame(interp_depths(SM_sim_mat, depths_SM), vars_SM)

  paths = [[:hydraulic, :kv, :kv, 2], [:hydraulic, :profile, :b, 2]]
  theta0 = parameters(model; paths).value
  gof = goodness_soilwater_θ1(theta0, model, forcing, dates;
    paths, θ1_obs=Float64.(θ1_obs), SM_obs, depths_SM=Float64.(depths_SM))

  @test nrow(gof) == length(depths_SM)
  @test all(isfinite, gof.KGE)
end
