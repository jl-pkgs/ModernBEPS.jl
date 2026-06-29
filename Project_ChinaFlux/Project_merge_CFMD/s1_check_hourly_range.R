#!/usr/bin/env Rscript
# hourly 值域与异常值核对 — 输出到 OUTPUT/hourly_value_check.csv
# 数据: Forcing_Hourly_Met_sp31_v20260614.csv
pacman::p_load(data.table, Ipaper, tibble)

f <- "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv"
d <- fread(f)

vars <- c("Ta_canopy", "RH_canopy", "WS_canopy", "Rs", "Rln_in", "Prcp")
# 合理值域（参照 Report §4.2 与 BEPS 物理范围）
# Rs:      0 ~ 1400 W/m²（地表太阳常数上限 ≈ 1361，地面 ~ 0~1100）
# Rln_in: 50 ~ 500 W/m²（夜间下行长波典型 200~400；负值/小正值为异常）
# 合理值域 + 异常阈值（异常指明显超出物理范围，如 Ta=-103.8, Uz=357）

# styler: off
bounds <- tribble(
  ~var,          ~lo,  ~hi, ~LOWER, ~UPPER,
  "Ta_canopy",   -40,   50,    -60,     60,
  "RH_canopy",     0,  100,    -10,    110,
  "WS_canopy",     0,   20,    -10,     30,
  "Rs",            0, 1400,    -10,   1500,
  "Rln_in",       50,  500,   -200,    800,
  "Prcp",          0,  100,    -10,    200
)
# styler: on

lo <- setNames(bounds$lo, bounds$var)
hi <- setNames(bounds$hi, bounds$var)
LOWER <- setNames(bounds$LOWER, bounds$var)
UPPER <- setNames(bounds$UPPER, bounds$var)

per_site <- d[,
  {
    out <- list(n = .N, nyear = length(unique(format(time, "%Y"))))
    for (v in vars) {
      x <- as.numeric(get(v))
      out[[paste0(v, "_NA_pct")]] <- round(100 * mean(!is.finite(x)), 2)
      x <- x[is.finite(x)]
      n <- length(x)
      out[[paste0(v, "_min")]] <- if (n) min(x) else NA_real_
      out[[paste0(v, "_max")]] <- if (n) max(x) else NA_real_
      out[[paste0(v, "_mean")]] <- if (n) mean(x) else NA_real_
      out[[paste0(v, "_p99")]] <- if (n) quantile(x, 0.99, names = FALSE) else NA_real_
      n_outlier <- sum(x < LOWER[v] | x > UPPER[v])
      out[[paste0(v, "_outlier_n")]] <- n_outlier
      out[[paste0(v, "_outlier_pct")]] <- round(100 * n_outlier / .N, 3)
      n_outofrange <- sum(x < lo[v] | x > hi[v])
      out[[paste0(v, "_outofrange_pct")]] <- round(100 * n_outofrange / .N, 3)
    }
    out
  },
  by = site
]

setorder(per_site, site)
out_file <- "./Project_ChinaFlux/Project_merge_CFMD/hourly_value_check.csv"
fwrite(per_site, out_file)

# 摘要：对每个变量，按"百分比列 > 阈值"筛出异常站点并打印
# 复用入口：pct_suffix / threshold / header_label / row_label 决定四种配置
report_anomalies <- function(per_site, vars, pct_suffix, threshold,
                             header_label, row_label) {
  cat(sprintf(
    "\n=== 摘要：%s (%s > %.2f%%) ===\n",
    header_label, sub("_pct$", "", pct_suffix), threshold
  ))

  for (v in vars) {
    pct_col <- paste0(v, pct_suffix)
    bad <- per_site[, .(site, pct = get(pct_col))][pct > threshold]
    if (nrow(bad) > 0) {
      ok(sprintf("\n[%s] %s的站点（%d 站）:\n", v, row_label, nrow(bad)))
      print(as_tibble(bad))
    }
  }
}

report_anomalies(per_site, vars, "_outofrange_pct", 0.5, "超合理值域", "超合理范围")
report_anomalies(per_site, vars, "_outlier_pct", 0.05, "异常值", "异常")
