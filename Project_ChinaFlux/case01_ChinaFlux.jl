## 1. 读取驱动数据
using BEPS, RTableTools, DataFrames, Dates, ModelParams, Ipaper
using JLD2
FT = Float64

indir = "z:/GitHub/jl-pkgs/ChinaFlux2026" |> path_mnt
st_full = fread("$indir/data/st_flux12.csv")


f = "$indir/data-raw/BEPS/Forcing_Met_Hourly_BEPS_Forest_sp12_hourly_v20260507.csv"
@time FORCING = fread("$indir/data-raw/BEPS/Forcing_Met_Hourly_BEPS_Forest_sp12_hourly_v20260507.csv")
replace_missing!(FORCING)
SITES_bad = ["MF_乔灌混交林_燕山"]
SITES = setdiff(unique(FORCING.site), SITES_bad)


##
function LoadData(SITE)
  # f = "$indir/data-raw/BEPS/Fluxes/$(SITE)_FluxLAISoil_daily_v20260510.csv" |> path_mnt
  f = "$indir/data-raw/Daily/BEPS/$(SITE)_FluxLAISoil_Daily_v20260511.csv" |> path_mnt

  FluxALL = fread(f)
  replace_missing!(FluxALL)
  rename!(FluxALL, :LAI_glass_G005 => :lai, :GPP => :GPP_obs, :ET => :ET_obs, :Hs => :Hs_obs)
  normalize_flux_obs!(FluxALL)
  (; lai) = FluxALL
  ntime2 = length(lai) * 24

  d_forcing = FORCING[FORCING.site.==SITE, 2:end]
  rename!(d_forcing, [:Ta_canopy => :Tair, :RH_canopy => :RH, :WS_canopy => :Uz])
  if (nrow(d_forcing) - ntime2) > 0
    @warn "驱动数据时间步长与率定数据不匹配，已截断驱动数据以匹配率定数据长度"
    d_forcing = d_forcing[1:ntime2, :]
  end
  clean_stats = sanitize_forcing!(d_forcing)
  @info "Forcing quality control" clean_stats
  (; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing
  ntime = length(Tair)
  forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)
  dates_local = parse_time.(d_forcing.time)

  dates_local, forcing, lai, FluxALL
end


function RunModel(SITE; maxn=1000, outdir="Project_ChinaFlux/OUTPUT",
  goal=:NSE, goal_multiplier=-1)

  fout = "$outdir/BEPS_$(SITE).jld2"
  printstyled("[site]: $SITE\n", color=:blue, bold=true, underline=true)
  # isfile(fout) && return

  # 率定所需数据
  # data-raw/Daily/BEPS/DBF_天然栎林_宝天曼_FluxLAISoil_Daily_v20260511.csv
  dates_local, forcing, lai, FluxALL = LoadData(SITE)
  dates_UTC = dates_local .- Hour(8) # [local] -> [UTC]

  ## 2. 初始化模型参数和状态变量
  st = st_full[findfirst(st_full.site .== SITE), :]
  (; lon, lat, VegType, SoilType, z_Uz, z_overstory) = st
  depths_SM = map(x -> parse(Int, x), split(st.z_SM, ",")) ./ 100 # 
  depths_TS = map(x -> parse(Int, x), split(st.z_TS, ",")) ./ 100 # 

  model = ParamBEPS(VegType, SoilType)
  model.veg.z_wind = z_Uz
  model.veg.z_canopy_o = z_overstory

  state = InitState0(model, forcing)
  @time df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates_UTC;
    ps=model, state, lon, lat)
  @time gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_local, FluxALL;
    depths_SM, depths_TS)
  gof
  # fwrite(cbind(; time=dates, df_ET = df_ET .* 3600), "./Project_ChinaFlux/df_ET.csv")
  
  ## 3. 参数优化模块
  opts = [
    (; path=[:r_drainage], name=:r_drainage, value=model.r_drainage),
    (; path=[:veg, :Ω], name=:Ω, value=model.veg.Ω),
    (; path=[:veg, :g1_w], name=:g1_w, value=model.veg.g1_w),
    (; path=[:veg, :g0_w], name=:g0_w, value=model.veg.g0_w),
    (; path=[:veg, :VCmax25], name=:VCmax25, value=model.veg.VCmax25),
  ] |> DataFrame
  paths = opts.path

  kw_loss = (; lon, lat, depths_SM, depths_TS, FluxDay=FluxALL,
    goal, goal_multiplier)

  @time theta_opt = optim(model, forcing, lai, dates_UTC; paths, maxn, kw_loss...)
  opts.theta_opt = theta_opt
  opts

  gof_opt, data_sim, data_obs = goodness(theta_opt, model, forcing, lai, dates_UTC; paths, kw_loss...)
  gof_opt
  jldsave(fout; gof_opt, gof, theta_opt, data_sim, data_obs)
end


for SITE in SITES[1:1]
  try
    RunModel(SITE; maxn=1000, outdir="Project_ChinaFlux/OUTPUT/NSE", goal=:NSE, goal_multiplier=-1)
  catch ex
    @error "Error processing site $SITE: $ex"
  end
end
# SITE = SITES[3]

# include("main_vis.jl")
# plot_result(data_sim, data_obs)
