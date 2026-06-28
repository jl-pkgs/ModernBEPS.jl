## 水量平衡诊断：逐站年降雨 P、年蒸发 ET(sim/obs)，验证降雨单位是否合理、找出无降雨站点
## 用法: julia --project Project_ChinaFlux/Plan/check_water_balance.jl
using JLD2, DataFrames, RTableTools, Dates, Statistics, Printf, Ipaper

indir = "z:/GitHub/jl-pkgs/ChinaFlux2026" |> path_mnt
FORCING = fread("$indir/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv")
replace_missing!(FORCING)
dir = "Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE"
files = sort(filter(f -> startswith(basename(f), "BEPS_") && endswith(f, ".jld2"), readdir(dir; join=true)))

out = DataFrame()
for f in files
  site = replace(basename(f), "BEPS_" => "", ".jld2" => "")
  d = FORCING[FORCING.site.==site, :]
  nrow(d) == 0 && continue
  nyear = nrow(d) / 8760
  P_ann = sum(skipmissing(d.Prcp)) / nyear          # 年降雨 [mm]
  o = load(f)
  ds, dob = o["data_sim"], o["data_obs"]
  ETs = mean(filter(isfinite, ds.ET)) * 365          # 年 ET sim [mm]
  ETo = mean(filter(isfinite, dob.ET)) * 365          # 年 ET obs [mm]
  push!(out, (; site,
    nyear=round(nyear, digits=1),
    P_ann=round(P_ann),
    ET_sim=round(ETs), ET_obs=round(ETo),
    ETsim_div_P=round(ETs / P_ann, digits=2),
    ETobs_div_P=round(ETo / P_ann, digits=2)); cols=:union)
end
fwrite(out, "$dir/water_balance.csv")

@printf("%-30s %5s %7s %7s %7s %7s %7s\n", "SITE", "yr", "P", "ETsim", "ETobs", "ETs/P", "ETo/P")
println("-"^78)
for r in eachrow(out)
  @printf("%-30s %5.1f %7.0f %7.0f %7.0f %7.2f %7.2f\n",
    r.site, r.nyear, r.P_ann, r.ET_sim, r.ET_obs, r.ETsim_div_P, r.ETobs_div_P)
end
println("\n无降雨站点 (P_ann < 50mm):")
foreach(s -> println("  ", s), out.site[out.P_ann.<50])
println("\n写出: $dir/water_balance.csv")
