using BEPS, Test, Dates, DataFrames

ntime = 6
forcing = MetSeries{Float64}(; ntime)
forcing.Tair .= 15.0
forcing.Prcp .= 1.0
dates = collect(DateTime(2020, 1, 1):Hour(1):DateTime(2020, 1, 1, ntime - 1))

model = ParamBEPS(25, 8)
state = InitState0(model, forcing)

##
@testset "soilwater standalone + calibration helpers" begin
  ntime = 6
  forcing = MetSeries{Float64}(; ntime)
  forcing.Tair .= 15.0
  forcing.Prcp .= 1.0
  dates = collect(DateTime(2020, 1, 1):Hour(1):DateTime(2020, 1, 1, ntime - 1))

  model = ParamBEPS(25, 8)
  state = InitState0(model, forcing)

  df = simulate_SM(forcing, dates; ps=model, state)
  @test nrow(df) == ntime
  @test "inf" in names(df)
  @test "z_water" in names(df)
  @test "θ1" in names(df)

  paths = [[:r_drainage]]
  theta = [model.r_drainage]
  df_pred = predict_SM(theta, model, forcing, dates; paths)
  @test nrow(df_pred) == ntime
  @test "inf" in names(df_pred)

  depths_SM = Float64[0.05, 0.15, 0.30, 0.60, 1.0]
  vars_SM = map(_sm_varname, depths_SM)
  nlayer = Int(model.N)
  θ_cols = [Symbol(:θ, j) for j in 1:nlayer]
  SM_sim_mat = hcat([df_pred[!, c] for c in θ_cols]...)
  SM_sim = interp_depths(SM_sim_mat, depths_SM)
  SM_obs = DataFrame()
  for j in eachindex(vars_SM)
    SM_obs[!, vars_SM[j]] = SM_sim[:, j]
  end

  gof = goodness_SM(theta, model, forcing, dates; paths, SM_obs, depths_SM)
  @test nrow(gof) == length(depths_SM)
  @test "KGE" in names(gof)
end
