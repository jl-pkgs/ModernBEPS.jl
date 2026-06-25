## 1. 读取驱动数据
using BEPS, RTableTools, DataFrames, Dates, ModelParams, Ipaper
using JLD2
FT = Float64

indir = "z:/GitHub/jl-pkgs/ChinaFlux2026" |> path_mnt
st_full = fread("$indir/data/Metadata/ChinaFlux_Metadata.csv")  # 39 站元数据，含 31 站全部


f = "$indir/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv"
@time FORCING = fread("$indir/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv")

replace_missing!(FORCING)
SITES_bad = String[]                                   # 按运行结果再补充异常站
SITES = setdiff(unique(FORCING.site), SITES_bad)        # 驱动文件即 31 站的权威清单


# 鲁棒列处理：仅重命名存在的列；保证列为 Float64（缺测→NaN），整列不存在→全 NaN
rename_existing!(d, pairs) = rename!(d, filter(p -> first(p) in propertynames(d), pairs))
# 元数据为手工维护：含 "NA" 等会把整列读成 String，统一转数值（不可解析→missing）
asnum(x) = x isa Number ? Float64(x) : (x isa AbstractString ? something(tryparse(Float64, x), missing) : missing)
ascol!(d, col) = d[!, col] = (col in propertynames(d)) ?
                             Float64.(coalesce.(d[!, col], NaN)) : fill(NaN, nrow(d))


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

  d_forcing = FORCING[FORCING.site.==SITE, 2:end]
  rename!(d_forcing, [:Ta_canopy => :Tair, :RH_canopy => :RH, :WS_canopy => :Uz])

  # 驱动可能比观测更早开始/更长（数据本身如此）：对齐到观测首日，再裁剪到观测长度
  i_beg = findfirst(t -> Date(parse_time(t)) == FluxALL.date[1], d_forcing.time)
  isnothing(i_beg) && error("forcing 不含观测首日 $(FluxALL.date[1])")
  d_forcing = d_forcing[i_beg:min(i_beg + ntime2 - 1, end), :]
  clean_stats = sanitize_forcing!(d_forcing)
  @info "Forcing quality control" clean_stats
  (; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing
  ntime = length(Tair)
  forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)
  dates_local = parse_time.(d_forcing.time)

  dates_local, forcing, lai, FluxALL
end


function RunModel(SITE; maxn=1000, outdir="Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE", overwrite=false,
  goal=:NSE, goal_multiplier=-1, SolveSM_fn=SolveSM_Bonan)

  mkpath(outdir)
  fout = "$outdir/BEPS_$(SITE).jld2"
  (isfile(fout) && !overwrite) && return

  printstyled("[site]: $SITE\n", color=:blue, bold=true, underline=true)

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
  parse_depths(s) = ismissing(s) ? Float64[] : map(x -> parse(Int, strip(x)), split(s, ",")) ./ 100
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

  jldsave(fout; gof_opt, gof, theta_opt, data_sim, data_obs)
  gof_opt
end



errors = Tuple{String,String}[]
outdir = "Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE"
for i in eachindex(SITES)
  !isCurrentWorker(i) && continue

  SITE = SITES[i]
  try
    RunModel(SITE; maxn=1000, outdir, goal=:NSE, goal_multiplier=-1, SolveSM_fn=SolveSM_Bonan)
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

# include("main_vis.jl")
# plot_result(data_sim, data_obs)
