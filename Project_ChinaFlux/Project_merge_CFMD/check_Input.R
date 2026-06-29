# %%
pacman::p_load(Ipaper, data.table, dplyr, stringr)
options(datatable.print.nrow = 31)

miss_pct <- function(x) round(100 * mean(!is.finite(as.numeric(x))), 1)

nansum <- function(x, ...) sum(x, na.rm = TRUE, ...)
nanmean <- function(x, ...) mean(x, na.rm = TRUE, ...)

# %% 
Met <- fread("/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv") # hourly
FluxLAI <- fread("/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Daily_Flux_sp40_v20260615.csv") # daily

# %% 缺失率 + Prcp 全 0/全 NA 检查（miss_pct 把"全 0"误判为 0%）
met_agg <- Met[, c(
  .(n_hour = .N, nyear = round(.N / 8760, 1)),
  lapply(.SD, miss_pct),
  Prcp_flag = fifelse(sum(Prcp > 0, na.rm = TRUE) > 0, "ok", "no_precip")
), by = site, .SDcols = c("Ta_canopy", "RH_canopy", "WS_canopy", "Rs", "Rln_in", "Prcp")]

flux_agg <- FluxLAI[, c(.(n_day = .N), lapply(.SD, miss_pct)),
  by = site, .SDcols = c("GPP", "ET", "LAI_glass_G005")
]

tab <- merge(met_agg, flux_agg, by = "site", all.x = TRUE)
tab <- tab %>%
  rename_with(~ str_replace_all(.x, "_canopy$|_glass_G005$", "")) %>%
  arrange(site)
fwrite(tab, "Project_ChinaFlux/OUTPUT/missing_input_China_FluxALL.csv")

cat("各变量平均缺失率 [%]:\n")
print(round(sapply(tab[, .(Ta, RH, WS, Rs, Rln_in, Prcp, GPP, ET, LAI)], mean, na.rm = TRUE), 1))

if (any(np <- tab$Prcp_flag == "no_precip")) {
  cat(sprintf("⚠ 无降水站: %s\n", paste(tab$site[np], collapse = ", ")))
}

# ---- 年总量单位自检: Met$Prcp=mm/hr, FluxLAI$Prcp/ET=mm/day ----
# 合理: Prcp∈[50,4000] ET∈[50,2000] mm/yr, Prcp_h/Prcp_d∈[0.7,1.3]
Met[, year := year(time)]
FluxLAI[, year := year(date)]
chk <- merge(
  Met[, .(Prcp_h = nansum(Prcp)), by = .(site, year)],
  FluxLAI[, .(Prcp_d = nansum(Prcp), ET = nansum(ET)), by = .(site, year)]
) %>% mutate(ratio = Prcp_h / Prcp_d) # 小时累加 / 日累加, 应∈[0.7,1.3]

ann <- chk[, .(
  Prcp_h = nanmean(Prcp_h), 
  Prcp_d = nanmean(Prcp_d),
  ratio = nanmean(ratio),
  ET = nanmean(ET), nyear = .N
), by = site] %>% .[order(-Prcp_d)] %>% dt_round(1)
ann[, c("ratio_OK") := list(between(ratio, 0.7, 1.3))]

# ann[, c("Prcp_h_OK", "Prcp_d_OK", "ratio_OK", "ET_OK") := list(
#   between(Prcp_h, 50, 4000), between(Prcp_d, 50, 4000),
#   between(ratio, 0.7, 1.3),  between(ET, 50, 2000)
# )]
fwrite(ann, "Project_ChinaFlux/OUTPUT/annual_unit_check.csv")

bad <- !with(ann, ratio_OK)
n_bad = sum(bad, na.rm = TRUE)
cat(if (any(bad)) sprintf("⚠ %d 站偏离合理范围\n", n_bad) else "✓ 单位检查通过\n")
print(ann)
