# %% 
pacman::p_load(
  Ipaper, data.table, dplyr, lubridate
)

# %% 
Met <- fread("/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv")
FluxLAI <- fread("/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Daily_Flux_sp40_v20260615.csv")

SITES = Met$site %>% unique_sort()
# SITE <- SITES[8]

# %% 数据缺失表：核对哪些站点/变量需要 ERA5-Land 补充
## 缺失比例 [%]：NA / NaN / Inf 均计为缺失；列不存在 -> 100%
miss_pct <- function(x) round(100 * mean(!is.finite(suppressWarnings(as.numeric(x)))), 1)

## 气象驱动：全部变量
met_vars  <- c("Ta_canopy", "RH_canopy", "WS_canopy", "Rs", "Rln_in", "Prcp")
met_miss  <- Met[, c(.(n_hour = .N, nyear = round(.N / 8760, 1)), lapply(.SD, miss_pct)), by = site, .SDcols = met_vars]

## FluxLAI：GPP, ET, LAI
flux_vars <- c("GPP", "ET", "LAI_glass_G005")
flux_miss <- FluxLAI[, c(.(n_day = .N), lapply(.SD, miss_pct)), by = site, .SDcols = flux_vars]

## 合并（以 31 个气象站为准），统一列名
tab <- merge(met_miss, flux_miss, by = "site", all.x = TRUE)
setnames(tab,
  c("Ta_canopy", "RH_canopy", "WS_canopy", "LAI_glass_G005"),
  c("Ta", "RH", "WS", "LAI"))
setcolorder(tab, c("site", "n_hour", "n_day", "nyear",
  "Ta", "RH", "WS", "Rs", "Rln_in", "Prcp", "GPP", "ET", "LAI"))
setorder(tab, site)

fout <- "Project_ChinaFlux/OUTPUT/missing_input_China_FluxALL.csv"
dir.create(dirname(fout), recursive = TRUE, showWarnings = FALSE)
fwrite(tab, fout)
cat(sprintf("已写出: %s （%d 站）\n", fout, nrow(tab)))
print(tab)

## 各变量平均缺失率，一眼看出最该补的变量
cat("\n各变量平均缺失率 [%]:\n")
vars <- c("Ta", "RH", "WS", "Rs", "Rln_in", "Prcp", "GPP", "ET", "LAI")
print(round(sapply(tab[, ..vars], mean, na.rm = TRUE), 1))
# %% 年总量单位自检（逐站点）
## 模型期望: Met$Prcp 单位 mm/hr（SoilWater.jl:43 显式 /3600/1000 换算到 m/s）
##           FluxLAI$Prcp、ET 单位 mm/day
## 物理合理区间 [mm/yr]: Prcp∈[50,4000]、ET∈[50,2000]、Prcp_h/Prcp_d∈[0.7,1.3]
## （若 Met$Prcp 误为 mm/day 重复 24h，ratio≈24；若误为 mm/半小时，ratio≈2）

Met[, year := year(time)]
FluxLAI[, year := year(date)]

prcp_h <- Met[, .(Prcp_h = sum(Prcp, na.rm = TRUE), nh = .N), by = .(site, year)]
prcp_d <- FluxLAI[, .(Prcp_d = sum(Prcp, na.rm = TRUE),
                      ET     = sum(ET,   na.rm = TRUE), nd = .N),
                  by = .(site, year)]
chk <- merge(prcp_h, prcp_d, by = c("site", "year"))[, ratio := Prcp_h / Prcp_d]

ann <- chk[, .(Prcp_h = mean(Prcp_h), Prcp_d = mean(Prcp_d),
               ratio  = mean(ratio, na.rm = TRUE),
               ET     = mean(ET, na.rm = TRUE), nyear = .N),
           by = site][order(-Prcp_d)]
ann[, c("Prcp_h_OK", "Prcp_d_OK", "ratio_OK", "ET_OK") := list(
  between(Prcp_h, 50, 4000), between(Prcp_d, 50, 4000),
  between(ratio,  0.7, 1.3), between(ET,   50, 2000))]

fout2 <- "Project_ChinaFlux/OUTPUT/annual_unit_check.csv"
fwrite(ann, fout2)

cat(sprintf("\n=== 年总量单位自检 (已写出 %s, %d 站) ===\n", fout2, nrow(ann)))
cat("模型期望: Met$Prcp=mm/hr, FluxLAI$Prcp/ET=mm/day\n")
cat("合理区间: Prcp∈[50,4000] ET∈[50,2000] mm/yr;  Prcp_h/Prcp_d∈[0.7,1.3]\n")
print(ann, nrow = nrow(ann))

bad <- ann[!(Prcp_h_OK & Prcp_d_OK & ratio_OK & ET_OK)]
if (nrow(bad)) {
  cat(sprintf("\n⚠ %d 站偏离合理范围（多为时段覆盖/数据质量问题）:\n", nrow(bad)))
  print(bad[, .(site, Prcp_h, Prcp_d, ratio, ET, nyear)], nrow = nrow(bad))
} else {
  cat("\n✓ 全部站点年总量在合理范围，单位检查通过\n")
}
