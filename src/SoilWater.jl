"""
    simulate_SM(forcing, dates; ps, state, ETi_obs, Tsoil_obs, kstep, SolveSM_fn)

单独运行土壤水分模块（不调用完整BEPS），用于快速调试和参数率定。

# Arguments
- `forcing`    : `MetSeries`，只使用 `Prcp` 字段 [mm/hr]
- `dates`      : 时间戳向量（长度应与 `forcing.ntime` 一致）
- `ps`         : `ParamBEPS`，土壤水参数
- `state`      : `StateBEPS` 初始状态（函数内部会 `deepcopy`）
- `ETi_obs`    : 各层蒸散发 [m s-1]，`nlayer × ntime` 矩阵；`nothing` 时全为0
- `Tsoil_obs`  : 各层土壤温度 [°C]，`nlayer × ntime`；`nothing` 时默认5°C（不冻结）
- `kstep`      : 子时间步长 [s]
- `SolveSM_fn` : 土壤水求解器函数，默认 `SolveSM_BEPS`

# Returns
DataFrame，列：`date`, `z_water`, `inf`, `θ1`…`θn`
"""
function simulate_SM(forcing::MetSeries, dates::AbstractVector;
  ps::ParamBEPS,
  state::StateBEPS,
  ETi_obs::Union{Nothing,AbstractMatrix}=nothing,
  Tsoil_obs::Union{Nothing,AbstractMatrix}=nothing,
  kstep::Float64=360.0,
  SolveSM_fn=SolveSM_BEPS)

  state = deepcopy(state)
  ntime = forcing.ntime
  kloop = round(Int, 3600.0 / kstep)
  n = Int(state.n_layer)

  θ_out = zeros(n, ntime)
  inf_out = zeros(ntime)
  zw_out = zeros(ntime)

  state.Tsoil_c[1:n] .= 5.0 # simplified warm-soil assumption (>0°C) to avoid freezing effects
  state.ETi[1:n] .= 0.0

  for i = 1:ntime
    !isnothing(Tsoil_obs) && (state.Tsoil_c[1:n] .= @view Tsoil_obs[:, i])
    !isnothing(ETi_obs) && (state.ETi[1:n] .= @view ETi_obs[:, i])

    state.r_rain_g = forcing.Prcp[i] / 3600.0 / 1000.0 # [mm/hr] -> [m/s]
    for _ = 1:kloop
      UpdateSoilMoisture(state, ps, kstep; fix_sm=false, SolveSM_fn)
    end

    θ_out[:, i] .= state.θ[1:n]
    inf_out[i] = state.inf
    zw_out[i] = state.z_water
  end

  df = DataFrame(:date => dates, :z_water => zw_out, :inf => inf_out)
  for j in 1:n
    df[!, Symbol("θ$j")] = θ_out[j, :]
  end
  df
end


_sm_varname(depth_m) = Symbol("SM_$(Int(depth_m * 100))cm")

"""
    predict_SM(theta, model, forcing, dates; paths, ETi_obs, Tsoil_obs, SolveSM_fn)

更新 `paths` 指定的水力参数后运行 `simulate_SM`，用于土壤水参数率定正演。
"""
function predict_SM(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, dates::Vector{DateTime};
  paths, ETi_obs=nothing, Tsoil_obs=nothing,
  SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  model = deepcopy(model)
  BEPS.update!(model, paths, theta)

  state = InitState0(model, forcing)
  simulate_SM(forcing, dates; ps=model, state, ETi_obs, Tsoil_obs, SolveSM_fn)
end


"""
    goodness_SM(theta, model, forcing, dates; paths, SM_obs, depths_SM, ...)

土壤水分拟合优度，返回含 KGE/NSE 等指标的 DataFrame。

# Arguments
- `SM_obs`    : 观测土壤湿度 DataFrame，列名为 `SM_Xcm`（X 为深度cm整数）
- `depths_SM` : 观测深度 [m] 向量，与模拟层对应
"""
function goodness_SM(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, dates::Vector{DateTime};
  paths, SM_obs::DataFrame, depths_SM::Vector{FT},
  ETi_obs=nothing, Tsoil_obs=nothing,
  SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  df_sim = predict_SM(theta, model, forcing, dates;
    paths, ETi_obs, Tsoil_obs, SolveSM_fn)

  vars_SM = map(_sm_varname, depths_SM)
  n = length(depths_SM)
  nlayer = Int(model.N)
  θ_cols = [Symbol(:θ, j) for j in 1:nlayer]

  SM_sim_mat = Matrix(df_sim[!, θ_cols])
  SM_sim = interp_depths(SM_sim_mat, depths_SM)

  gof = map(1:n) do j
    obs_col = SM_obs[!, vars_SM[j]]
    sim_col = SM_sim[:, j]
    mask = map(x -> !ismissing(x) && !isnan(Float64(x)), obs_col)
    (; var=vars_SM[j], GOF(Float64.(obs_col[mask]), sim_col[mask])...)
  end |> DataFrame
  gof
end


"""
    optim_SM(model, forcing, dates; paths, SM_obs, depths_SM, maxn=200, ...)

仅针对土壤水力参数的自动率定，比全模型率定快约10倍。

# Arguments
- `paths`     : 待率定参数路径，例如 `[(:hydraulic, :K_sat), (:hydraulic, :b)]`
- `SM_obs`    : 观测土壤湿度 DataFrame，列名为 `SM_Xcm`
- `depths_SM` : 观测深度 [m] 向量
- `goal`      : 目标函数，`:KGE` 或 `:NSE`（默认 `:KGE`）
- `maxn`      : SCE-UA 最大迭代次数
"""
function optim_SM(model::ParamBEPS{FT}, forcing::MetSeries{FT},
  dates::Vector{DateTime};
  paths=nothing,
  SM_obs::DataFrame, depths_SM::Vector{FT},
  goal::Symbol=:KGE, maxn::Int=200,
  ETi_obs=nothing, Tsoil_obs=nothing,
  SolveSM_fn=SolveSM_BEPS) where {FT<:AbstractFloat}

  params_all = parameters(model;)
  isnothing(paths) && (paths = filter(x -> x.path[1] == :hydraulic, params_all))

  function _loss(theta)
    gof = goodness_SM(theta, model, forcing, dates;
      paths, SM_obs, depths_SM, ETi_obs, Tsoil_obs, SolveSM_fn)
    -mean(gof[!, goal])
  end

  params = parameters(model; paths)
  lb = FT.(first.(params.bound))
  ub = FT.(last.(params.bound))
  u0 = params.value
  theta, _, _ = sceua(_loss, u0, lb, ub; maxn, verbose=true, parallel=true)
  theta
end

export predict_SM, goodness_SM, optim_SM
export simulate_SM
export _sm_varname
