using JLD2, Printf

vers = ["NSE_CMFD_1h", "NSE_CMFD_1h_V2", "NSE_CMFD_1h_V3"]
base = "Project_ChinaFlux/OUTPUT/ALL/Bonan"

function read_gof(dir, site)
  f = joinpath(base, dir, "BEPS_$site.jld2")
  isfile(f) || return missing
  jldopen(f) do d
    g = d["gof_opt"]
    (GPP_NSE = g.GPP_NSE, ET_NSE = g.ET_NSE, Hs_NSE = g.Hs_NSE, FluxNSE = g.FluxNSE)
  end
end

files = readdir(joinpath(base, vers[end]))
sites = replace.(filter(f -> endswith(f, ".jld2"), files), "BEPS_" => "", ".jld2" => "")

io = open("Project_ChinaFlux/OUTPUT/ALL/Bonan/V3_comparison.txt", "w")
println(io, rpad("站点", 44), "  V1_GPP  V1_ET  V2_GPP  V2_ET  V3_GPP  V3_ET   ΔvsV2")
println(io, repeat("-", 120))
for s in sort(sites)
  g1 = read_gof(vers[1], s)
  g2 = read_gof(vers[2], s)
  g3 = read_gof(vers[3], s)
  ismissing(g3) && continue
  vn(n) = ismissing(n) ? NaN : round(Float64(n), digits=3)
  v1g, v1e = vn(g1.GPP_NSE), vn(g1.ET_NSE)
  v2g, v2e = vn(g2.GPP_NSE), vn(g2.ET_NSE)
  v3g, v3e = vn(g3.GPP_NSE), vn(g3.ET_NSE)
  refg = isnan(v2g) ? v1g : v2g; refe = isnan(v2e) ? v1e : v2e
  dg, de = round(v3g - refg, digits=3), round(v3e - refe, digits=3)
  sg = dg > 0 ? "↑" : dg < 0 ? "↓" : "="
  se = de > 0 ? "↑" : de < 0 ? "↓" : "="
  @printf(io, "%-44s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f  GPP%s%.3f ET%s%.3f\n",
    s, v1g, v1e, v2g, v2e, v3g, v3e, sg, abs(dg), se, abs(de))
end
close(io)
println("Done: $(joinpath(base, "V3_comparison.txt"))")
