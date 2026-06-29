using DataFrames, RTableTools, Printf, Statistics

function load_flux_nse(dir)
  f = joinpath(dir, "GOF_summary.csv")
  isfile(f) || return Dict()
  d = fread(f)
  d = d[d.var .∈ Ref(["GPP", "ET"]), :]
  tof(x) = x isa Number ? Float64(x) : something(tryparse(Float64, string(x)), NaN)
  vals = Dict{String, Tuple{Float64, Float64}}()
  for site in unique(d.site)
    g = d[(d.site .== site) .& (d.var .== "GPP"), :NSE_opt]
    e = d[(d.site .== site) .& (d.var .== "ET"), :NSE_opt]
    vals[site] = (isempty(g) ? NaN : tof(g[1]), isempty(e) ? NaN : tof(e[1]))
  end
  vals
end

mean_flux(g, e) = (x = filter(!isnan, [g, e]); isempty(x) ? NaN : mean(x))

base = "Project_ChinaFlux/OUTPUT/ALL/Bonan"
v2 = load_flux_nse(joinpath(base, "NSE_CMFD_1h_V2"))
v4 = load_flux_nse(joinpath(base, "NSE_CMFD_1h_V4"))

n_up, n_dn, n_eq = 0, 0, 0
@printf("%-28s %7s %7s %7s %7s  %s\n", "站点", "V2_GPP", "V2_ET", "V4_GPP", "V4_ET", "FluxΔ")
println(repeat("-", 85))
for site in sort(collect(keys(v4)))
  g2, e2 = v2[site]; g4, e4 = v4[site]
  f2 = mean_flux(g2, e2); f4 = mean_flux(g4, e4)
  df = isnan(f2) ? NaN : round(f4 - f2, digits = 3)
  if df > 0.005; global n_up += 1
  elseif df < -0.005; global n_dn += 1
  else global n_eq += 1
  end
  s = df > 0.005 ? "↑" : df < -0.005 ? "↓" : "="
  @printf("%-28s %7.3f %7.3f %7.3f %7.3f  %s%.3f\n", site, g2, e2, g4, e4, s, abs(df))
end
println(repeat("-", 85))
@printf("合计: ↑%d ↓%d =%d\n", n_up, n_dn, n_eq)
