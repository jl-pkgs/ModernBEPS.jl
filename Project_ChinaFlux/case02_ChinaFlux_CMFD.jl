## 1. 读取驱动数据
using BEPS, RTableTools, DataFrames, Dates, ModelParams, Ipaper
using JLD2
FT = Float64

indir = "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026"
st_full = fread("$indir/data/Metadata/ChinaFlux_Metadata.csv")  # 39 站元数据，含 31 站全部


f = "/mnt/z/China/CMFD_V2.0/OUTPUT/ChinaFlux_sp38_final_forcing_1h.csv"
@time FORCING = fread(f)
replace_missing!(FORCING)
rename_existing!(d, pairs) = rename!(d, filter(p -> first(p) in propertynames(d), pairs))
rename_existing!(FORCING, [
  :datetime => :time,
  :Temp => :Tair,
  :RHum => :RH,
  :Wind => :Uz,
  :SRad => :Rs,
  :LRad => :Rln_in,
  :Prec => :Prcp,
])

flux_files = filter(f -> endswith(f, "_Daily_FluxALL_v20260615.csv"),
  readdir("$indir/data/BEPS/Daily_FluxALL"))
SITES_obs = replace.(flux_files, "_Daily_FluxALL_v20260615.csv" => "")
SITES = sort(intersect(unique(FORCING.site), SITES_obs))
SITES_missing_obs = sort(setdiff(unique(FORCING.site), SITES_obs))
SITES_missing_forcing = sort(setdiff(SITES_obs, unique(FORCING.site)))
isempty(SITES_missing_obs) || @warn "CMFD forcing sites without matching FluxALL daily file" SITES_missing_obs
isempty(SITES_missing_forcing) || @warn "FluxALL daily sites without matching CMFD forcing" SITES_missing_forcing
if haskey(ENV, "BEPS_SITES")
  SITES_requested = strip.(split(ENV["BEPS_SITES"], ","))
  SITES = intersect(SITES, SITES_requested)
  isempty(SITES) && error("BEPS_SITES 与可运行站点无交集")
end
env_bool(name, default=false) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "y")
MAXN = parse(Int, get(ENV, "BEPS_MAXN", "1000"))
OVERWRITE = env_bool("BEPS_OVERWRITE")

# PRCP_SCALE = Dict(
#   "CRO_冬小麦夏玉米_固城" => 1 / 29.1,  # h_max 430mm/hr，年 15791mm → ~543mm
#   "CRO_水稻_盘锦" => 1 / 2.0,           # 小时为日累计的 2 倍
#   "CRO_水稻_句容" => 1 / 4.5,           # 小时为日累计的 4~5 倍
# )


# 鲁棒列处理：仅重命名存在的列；保证列为 Float64（缺测→NaN），整列不存在→全 NaN
# 元数据为手工维护：含 "NA" 等会把整列读成 String，统一转数值（不可解析→missing）
asnum(x) = x isa Number ? Float64(x) : (x isa AbstractString ? something(tryparse(Float64, x), missing) : missing)
ascol!(d, col) = d[!, col] = (col in propertynames(d)) ?
                             Float64.(coalesce.(d[!, col], NaN)) : fill(NaN, nrow(d))
cmfd_time(t::DateTime) = t
cmfd_time(t) = DateTime(first(String(t), 19), dateformat"yyyy-mm-ddTHH:MM:SS")
cmfd_local_time(t) = cmfd_time(t) + Hour(8)


##
function LoadData(SITE)
  f = "$indir/data/BEPS/Daily_FluxALL/$(SITE)_Daily_FluxALL_v20260615.csv" |> path_mnt

  FluxALL = fread(f)
  replace_missing!(FluxALL)
  # 部分站缺碳通量列（GPP/Hs），rename 仅作用于存在的列
  rename_existing!(FluxALL, [:LAI_glass_G005 => :lai, :GPP => :GPP_obs, :ET => :ET_obs, :Hs => :Hs_obs])
  foreach(c -> ascol!(FluxALL, c), [:lai, :GPP_obs, :ET_obs, :Hs_obs])  # 缺列填 NaN，GOF 自然返回 NaN
  normalize_flux_obs!(FluxALL)
  # 观测已为标准日尺度单位（GPP: gC m⁻² d⁻¹, ET: mm d⁻¹, Hs: W m⁻²），与 agg_daily 输出一致，无需换算
  # （单位见 data/BEPS/BEPS_Forcing_China_FluxALL.md §1.2）
  (; lai) = FluxALL
  ntime2 = length(lai) * 24

  d_forcing = FORCING[FORCING.site.==SITE, :]
  sort!(d_forcing, :time)

  # 驱动可能比观测更早开始/更长（数据本身如此）：对齐到观测首日，再裁剪到观测长度
  i_beg = findfirst(t -> Date(cmfd_local_time(t)) == FluxALL.date[1], d_forcing.time)
  isnothing(i_beg) && error("forcing 不含观测首日 $(FluxALL.date[1])")
  d_forcing = d_forcing[i_beg:min(i_beg + ntime2 - 1, end), :]
  # 降水异常已在上游数据修复，无需再按站缩放（保留以备回退）：
  # haskey(PRCP_SCALE, SITE) && (d_forcing.Prcp .*= PRCP_SCALE[SITE])  # 逐站降雨校正
  clean_stats = sanitize_forcing!(d_forcing)
  @info "Forcing quality control" clean_stats
  (; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing
  ntime = length(Tair)
  forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)
  dates_local = cmfd_local_time.(d_forcing.time)

  dates_local, forcing, lai, FluxALL
end


function RunModel(SITE; maxn=1000, outdir="Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE_CMFD_1h", overwrite=false,
  goal=:NSE, goal_multiplier=-1, SolveSM_fn=SolveSM_Bonan)

  mkpath(outdir)
  fout = "$outdir/BEPS_$(SITE).jld2"
  (isfile(fout) && !overwrite) && return

  printstyled("[site]: $SITE\n", color=:blue, bold=true, underline=true)
  t_beg = time()  # 记录单站运行时长（主要由 optim 决定），存入结果

  # 率定所需数据
  dates_local, forcing, lai, FluxALL = LoadData(SITE)
  dates_UTC = dates_local .- Hour(8) # [local] -> [UTC]

  ## 2. 初始化模型参数和状态变量
  i_st = findfirst(st_full.site .== SITE)
  isnothing(i_st) && error("missing metadata: $SITE")
  st = st_full[i_st, :]
  (; lon, lat, VegType, SoilType) = st
  z_Uz, z_overstory = asnum(st.z_Uz), asnum(st.z_overstory)  # 元数据手工维护，统一转数值
  SoilType = ismissing(SoilType) ? "loam" : SoilType    # 元数据缺失兜底

  # 元数据 z_SM/z_TS 可能缺失；缺失→不评价土壤
  parse_depths(s) = (ismissing(s) || isempty(strip(String(s)))) ?
                    Float64[] : map(x -> parse(Int, strip(x)), split(String(s), ",")) ./ 100
  depths_SM = parse_depths(st.z_SM)
  depths_TS = parse_depths(st.z_TS)

  # 深度过滤：各站日文件土壤列高度异构、命名可能不规范（如固城 z_SM=4 但列名为 SM_4cm_N）。
  # 仅保留标准列 SM_{d}cm / TS_{d}cm 存在「且含有限观测」的深度；无匹配则该站只评价 GPP/ET/Hs。
  function has_obs(prefix, d)
    c = Symbol("$(prefix)_$(Int(round(d * 100)))cm")
    c in propertynames(FluxALL) && any(isfinite, coalesce.(FluxALL[!, c], NaN))
  end
  depths_SM = filter(d -> has_obs("SM", d), depths_SM)
  depths_TS = filter(d -> has_obs("TS", d), depths_TS)

  model = ParamBEPS(VegType, SoilType)
  ismissing(z_Uz) || (model.veg.z_wind = z_Uz)          # z 高度缺失→保留模型默认
  ismissing(z_overstory) || (model.veg.z_canopy_o = z_overstory)

  state = InitState0(model, forcing)
  @time df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates_UTC;
    ps=model, state, lon, lat, SolveSM_fn)
  @time gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_local, FluxALL;
    depths_SM, depths_TS)

  ## 3. 参数优化（SCE-UA, maxn=1000）
  opts = [
    (; path=[:r_drainage], name=:r_drainage, value=model.r_drainage),
    (; path=[:veg, :Ω], name=:Ω, value=model.veg.Ω),
    (; path=[:veg, :g1_w], name=:g1_w, value=model.veg.g1_w),
    (; path=[:veg, :g0_w], name=:g0_w, value=model.veg.g0_w),
    (; path=[:veg, :VCmax25], name=:VCmax25, value=model.veg.VCmax25),
  ] |> DataFrame
  paths = opts.path

  kw_loss = (; lon, lat, depths_SM, depths_TS, FluxDay=FluxALL,
    goal, goal_multiplier, SolveSM_fn)

  @time theta_opt = optim(model, forcing, lai, dates_UTC; paths, maxn, kw_loss...)
  opts.theta_opt = theta_opt
  gof_opt, data_sim, data_obs = goodness(theta_opt, model, forcing, lai, dates_UTC; paths, kw_loss...)

  runtime = round(time() - t_beg, digits=1)  # [秒]
  jldsave(fout; gof_opt, gof, theta_opt, data_sim, data_obs, runtime)
  gof_opt
end



errors = Tuple{String,String}[]
outdir = "Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE_CMFD_1h"
for i in eachindex(SITES)
  !isCurrentWorker(i) && continue

  SITE = SITES[i]
  try
    RunModel(SITE; maxn=MAXN, outdir, goal=:NSE, goal_multiplier=-1,
      SolveSM_fn=SolveSM_Bonan, overwrite=OVERWRITE)
  catch ex
    msg = sprint(showerror, ex)
    @error "Error processing site $SITE: $msg"
    push!(errors, (SITE, msg))
  end
end

# 错误清单：便于「全部跑完后」按站定位问题
if !isempty(errors)
  @warn "以下站点运行失败" errors
  fwrite(DataFrame(site=first.(errors), error=last.(errors)), "$outdir/_errors.csv")
end


## 差站重跑（overwrite=true）
# 分级标记为「差」的站点（见 Plan/Report_China_FluxALL.typ §3）：使用 CMFD 强制重算。
SITES_poor = [
  "CRO_冬小麦夏玉米_固城",
  "CRO_水稻_盘锦",
  "CRO_水稻_长岭",
  "CRO_水稻_句容",
  "GRA_高寒草甸_若尔盖",
  "GRA_人工垂穗披碱草_三江源",
]
SITES_poor_run = intersect(SITES_poor, SITES)
for i in eachindex(SITES_poor_run)
  !isCurrentWorker(i) && continue
  SITE = SITES_poor_run[i]
  try
    RunModel(SITE; maxn=MAXN, outdir, goal=:NSE, goal_multiplier=-1,
      SolveSM_fn=SolveSM_Bonan, overwrite=true)
  catch ex
    @error "Error processing site $SITE: $(sprint(showerror, ex))"
  end
end

# include("main_vis.jl")
# plot_result(data_sim, data_obs)
