export InitState0, predict, goodness, loss, optim, BEPS_GOF

using ModelParams: sceua, of_KGE, of_NSE, GOF
import DataFrames: AbstractDataFrame

include("check_forcing.jl")

function InitState0(model::ParamBEPS{FT}, forcing::MetSeries{FT}) where {FT<:AbstractFloat}
  Ta = forcing.Tair[1]
  Tsoil0 = Ta
  θ0 = model.hydraulic.θ_sat[1] * 0.8  # 初始化为田间持水量（避免低于凋萎含水量）
  z_snow0 = 0.0
  BEPS._init_state(model, Tsoil0, Ta, θ0, z_snow0) # state
end


function predict(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, lai::Vector{FT}, dates_UTC::Vector{DateTime};
  paths, lon::FT, lat::FT, SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  model = deepcopy(model)
  params = parameters(model; paths)
  theta_prev = params.value
  BEPS.update!(model, paths, theta)

  state = InitState0(model, forcing)
  df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates_UTC;
    ps=model, state, lon, lat, SolveSM_fn)

  BEPS.update!(model, paths, theta_prev) # 恢复原参数值
  df_fluxes, df_ET, states, caches
end

function goodness(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, lai::Vector{FT}, dates_UTC::Vector{DateTime};
  paths, lon::FT, lat::FT, depths_SM::Vector{FT}, depths_TS::Vector{FT}, FluxDay::DataFrame,
  SolveSM_fn=SolveSM_BEPS, ignored...) where {FT<:AbstractFloat}

  df_fluxes, df_ET, states, caches = predict(theta, model, forcing, lai, dates_UTC;
    paths, lon, lat, SolveSM_fn)

  dates_local = dates_UTC .+ Hour(8) # [UTC] -> [local]
  gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_local, FluxDay;
    depths_SM, depths_TS)
  return gof, data_sim, data_obs
end

# dates: UTC时间，需转换为当地时间
function loss(theta::Vector{FT}, model::ParamBEPS{FT}, 
  forcing::MetSeries{FT}, lai::Vector{FT}, dates_UTC::Vector{DateTime}; 
  goal = :NSE, goal_multiplier = -1, 
  kw_loss...) where {FT<:AbstractFloat}
  gof, data_sim, data_obs = goodness(theta, model, forcing, lai, dates_UTC; kw_loss...)
  mean(gof.Flux[1:2, goal]) * goal_multiplier
end


# `dates`: 注意是UTC时间
"""
- `kw`: 
"""
function optim(model, forcing, lai, dates_UTC; paths, maxn=200, kw_loss...)
  _loss(theta) = loss(theta, model, forcing, lai, dates_UTC; paths, kw_loss...)

  params = parameters(model; paths)
  lb = map(x -> FT(x[1]), params.bound)
  ub = map(x -> FT(x[2]), params.bound)

  u0 = params.value
  theta, feval, exitflag = sceua(_loss, u0, lb, ub; maxn, verbose=true, parallel=true)
  theta
end



# GLOBAL: depths_SM, depths_TS
function BEPS_GOF(df_fluxes, states, dates_hour, FluxALL; depths_SM, depths_TS)
  ## 模拟
  FluxALL = copy(FluxALL)
  normalize_flux_obs!(FluxALL)
  obs_dates = _daily_obs_dates(FluxALL)
  (; lai, GPP_obs, ET_obs, Hs_obs) = FluxALL
  vars_TS = map(i -> Symbol("TS_$(Int(depths_TS[i] * 100))cm"), eachindex(depths_TS))
  vars_SM = map(i -> Symbol("SM_$(Int(depths_SM[i] * 100))cm"), eachindex(depths_SM))
  TS_obs = FluxALL[:, vars_TS] |> Matrix
  SM_obs = FluxALL[:, vars_SM] |> Matrix

  ## 观测
  (dates, GPP, ET, Hs) = agg_daily(df_fluxes, dates_hour)

  ## SM GOF
  SM_sim = states.vectors.θ
  TS_sim = states.vectors.Tsoil_c

  SM_day = agg_daily(SM_sim, dates_hour)
  TS_day = agg_daily(TS_sim, dates_hour)

  ## SM & TS
  SM = interp_depths(SM_day, depths_SM)
  TS = interp_depths(TS_day, depths_TS)

  sim_dates = dates
  obs_dates = isnothing(obs_dates) ? sim_dates : obs_dates
  data_sim = hcat(
    DataFrame(; date=sim_dates, GPP, ET, Hs),
    DataFrame(SM, vars_SM),
    DataFrame(TS, vars_TS))
  data_obs = hcat(
    DataFrame(; date=obs_dates, GPP=GPP_obs, ET=ET_obs, Hs=Hs_obs),
    DataFrame(SM_obs, vars_SM),
    DataFrame(TS_obs, vars_TS))

  if length(sim_dates) != length(obs_dates) || sim_dates != obs_dates
    data_sim, data_obs = _align_daily_data(data_sim, data_obs)
  end

  gof_SM = map(i -> (; var=vars_SM[i], GOF(data_obs[!, vars_SM[i]], data_sim[!, vars_SM[i]])...),
    eachindex(depths_SM)) |> DataFrame
  gof_TS = map(i -> (; var=vars_TS[i], GOF(data_obs[!, vars_TS[i]], data_sim[!, vars_TS[i]])...),
    eachindex(depths_TS)) |> DataFrame

  gof_Flux = DataFrame([
    (; var="GPP", GOF(data_obs.GPP, data_sim.GPP)...),
    (; var="ET", GOF(data_obs.ET, data_sim.ET)...),
    (; var="Hs", GOF(data_obs.Hs, data_sim.Hs)...)
  ])

  gof = (; Flux=gof_Flux, SM=gof_SM, TS=gof_TS)
  gof, data_sim, data_obs
end


"""
    predict_soilwater_θ1(theta, model, forcing, dates; paths, θ1_obs, ETi_obs, Tsoil_obs, SolveSM_fn)

以第1层观测 θ 为上边界驱动，更新 `paths` 参数后运行 `simulate_soilwater`。
"""
function predict_soilwater_θ1(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, dates::Vector{DateTime};
  paths,
  θ1_obs::AbstractVector{FT},
  ETi_obs=nothing, Tsoil_obs=nothing,
  SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  model = deepcopy(model)
  theta_prev = parameters(model; paths).value
  BEPS.update!(model, paths, theta)

  state = InitState0(model, forcing)
  df = simulate_soilwater(forcing, dates;
    ps=model, state, ETi_obs, Tsoil_obs, θ1_obs, SolveSM_fn)

  BEPS.update!(model, paths, theta_prev)
  df
end

"""
    goodness_soilwater_θ1(theta, model, forcing, dates; paths, θ1_obs, SM_obs, depths_SM, ...)

以第1层观测 θ 为上边界，计算第 2~n 层的拟合优度（KGE/NSE 等）。
"""
function goodness_soilwater_θ1(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, dates::Vector{DateTime};
  paths,
  θ1_obs::AbstractVector{FT},
  SM_obs::DataFrame,
  depths_SM::Vector{FT},
  ETi_obs=nothing, Tsoil_obs=nothing,
  SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  df_sim = predict_soilwater_θ1(theta, model, forcing, dates;
    paths, θ1_obs, ETi_obs, Tsoil_obs, SolveSM_fn)

  vars_SM = map(i -> Symbol("SM_$(Int(depths_SM[i] * 100))cm"), eachindex(depths_SM))
  n = length(depths_SM)

  θ_cols = filter(name -> startswith(String(name), "θ"), propertynames(df_sim))
  sort!(θ_cols; by=name -> parse(Int, replace(String(name), "θ" => "")))
  SM_sim_mat = hcat([df_sim[!, col] for col in θ_cols]...) # ntime × nlayer
  SM_sim = interp_depths(FT.(SM_sim_mat), depths_SM) # ntime × n

  gof = map(1:n) do j
    obs_col = SM_obs[!, vars_SM[j]]
    sim_col = SM_sim[:, j]
    mask = map(obs_col) do x
      !ismissing(x) && !isnan(Float64(x))
    end
    (; var=vars_SM[j], GOF(Float64.(obs_col[mask]), sim_col[mask])...)
  end |> DataFrame

  gof
end

"""
    optim_soilwater_θ1(model, forcing, dates; paths, θ1_obs, SM_obs, depths_SM, ...)

以第1层观测 θ 为上边界驱动，率定深层土壤水力参数。
"""
function optim_soilwater_θ1(model::ParamBEPS{FT}, forcing::MetSeries{FT},
  dates::Vector{DateTime};
  paths,
  θ1_obs::AbstractVector{FT},
  SM_obs::DataFrame,
  depths_SM::Vector{FT},
  goal::Symbol=:KGE, maxn::Int=200,
  ETi_obs=nothing, Tsoil_obs=nothing,
  SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  function _loss(theta)
    gof = goodness_soilwater_θ1(theta, model, forcing, dates;
      paths, θ1_obs, SM_obs, depths_SM, ETi_obs, Tsoil_obs, SolveSM_fn)
    -mean(gof[!, goal])
  end

  params = parameters(model; paths)
  lb = map(x -> FT(x[1]), params.bound)
  ub = map(x -> FT(x[2]), params.bound)
  u0 = params.value

  theta, feval, exitflag = sceua(_loss, u0, lb, ub; maxn, verbose=true, parallel=true)
  theta
end

export predict_soilwater_θ1, goodness_soilwater_θ1, optim_soilwater_θ1
