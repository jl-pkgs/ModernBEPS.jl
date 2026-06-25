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
met_miss  <- Met[, c(.(n_hour = .N), lapply(.SD, miss_pct)), by = site, .SDcols = met_vars]

## FluxLAI：GPP, ET, LAI
flux_vars <- c("GPP", "ET", "LAI_glass_G005")
flux_miss <- FluxLAI[, c(.(n_day = .N), lapply(.SD, miss_pct)), by = site, .SDcols = flux_vars]

## 合并（以 31 个气象站为准），统一列名
tab <- merge(met_miss, flux_miss, by = "site", all.x = TRUE)
setnames(tab,
  c("Ta_canopy", "RH_canopy", "WS_canopy", "LAI_glass_G005"),
  c("Ta", "RH", "WS", "LAI"))
setcolorder(tab, c("site", "n_hour", "n_day",
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
