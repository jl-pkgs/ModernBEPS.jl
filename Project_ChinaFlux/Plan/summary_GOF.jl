## 汇总已完成站点的拟合优度（优化前 gof / 优化后 gof_opt）
## 用法: julia --project Project_ChinaFlux/summary_GOF.jl
using JLD2, DataFrames, RTableTools, Printf

dir = "Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE"
files = sort(filter(f -> startswith(basename(f), "BEPS_") && endswith(f, ".jld2"),
  readdir(dir; join=true)))

# 取某指标 DataFrame 中某变量的一行（R² 列名兼容 :R2 / :R²）
getval(r, k) = hasproperty(r, k) ? getproperty(r, k) : NaN
# GPP/ET 缺测时 GOF 返回 -999 哨兵（有效样本 ≤ 2）或 NaN（序列退化）；汇总成表时统一转 NA（missing），避免误读为有效精度
na(x) = (ismissing(x) || (x isa Real && (isnan(x) || x == -999.0))) ? missing : x
function gofrow(g, var)
  i = findfirst(==(var), g.var)
  isnothing(i) && return (KGE=missing, NSE=missing, R2=NaN, RMSE=NaN, bias=NaN)
  r = g[i, :]
  (; KGE=na(getval(r, :KGE)), NSE=na(getval(r, :NSE)),
    R2=getval(r, hasproperty(r, :R2) ? :R2 : Symbol("R²")),
    RMSE=getval(r, :RMSE), bias=getval(r, :bias))
end

long = DataFrame()
runtimes = DataFrame()
for f in files
  site = replace(basename(f), "BEPS_" => "", ".jld2" => "")
  jldopen(f) do o
    (haskey(o, "gof") && haskey(o, "gof_opt")) || return
    push!(runtimes, (; site, runtime_s=haskey(o, "runtime") ? o["runtime"] : missing); cols=:union)
    for (stage, gof) in (("0", o["gof"]), ("opt", o["gof_opt"]))
      for grp in (:Flux, :SM, :TS)
        g = getfield(gof, grp)
        for r in eachrow(g)
          m = gofrow(g, r.var)
          push!(long, (; site, group=string(grp), var=r.var, stage,
              KGE=m.KGE, NSE=m.NSE, R2=m.R2, RMSE=m.RMSE, bias=m.bias); cols=:union)
        end
      end
    end
  end
end

fwrite(long, "$dir/GOF_summary_long.csv"; missingstring="NA")  # 缺测 NSE/KGE 写为 NA

# Flux 宽表：每站每变量一行，优化前后 KGE/NSE 并列
flux = long[long.group.=="Flux", :]
wide = DataFrame()
for sub in groupby(flux, [:site, :var])
  g0 = sub[sub.stage.=="0", :]
  g1 = sub[sub.stage.=="opt", :]
  push!(wide, (; site=sub.site[1], var=sub.var[1],
      KGE0=round(g0.KGE[1], digits=3), KGE_opt=round(g1.KGE[1], digits=3),
      NSE0=round(g0.NSE[1], digits=3), NSE_opt=round(g1.NSE[1], digits=3)); cols=:union)
end
if "runtime_s" in names(wide) || nrow(runtimes) > 0
  wide = leftjoin(wide, runtimes, on=:site)
end
fwrite(wide, "$dir/GOF_summary.csv"; missingstring="NA")  # 缺测 NSE/KGE 写为 NA
fwrite(runtimes, "$dir/runtime.csv")

println("已完成站点: ", length(unique(long.site)), " 站\n")
@printf("%-32s %-4s | %6s %6s | %6s %6s\n", "SITE", "var", "KGE0", "KGE*", "NSE0", "NSE*")
println("-"^66)
fmt(x) = (ismissing(x) || isnan(x)) ? "    NA" : @sprintf("%6.2f", x)
for r in eachrow(wide)
  @printf("%-32s %-4s | %6s %6s | %6s %6s\n",
    r.site, r.var, fmt(r.KGE0), fmt(r.KGE_opt), fmt(r.NSE0), fmt(r.NSE_opt))
end
println("\n写出: $dir/GOF_summary.csv (Flux 宽表), GOF_summary_long.csv (全指标)")
