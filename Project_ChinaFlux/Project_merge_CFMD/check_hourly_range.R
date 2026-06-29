#!/usr/bin/env Rscript
# hourly 值域与哨兵值核对 — 输出到 OUTPUT/hourly_value_check.csv
# 数据: Forcing_Hourly_Met_sp31_v20260614.csv
# 列名: site, time, Ta_canopy, RH_canopy, WS_canopy, Rs, Rln_in, Prcp

suppressPackageStartupMessages(library(data.table))

f <- "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv"
cat("Reading", f, "\n")
d <- fread(f)
cat("Loaded", nrow(d), "rows ×", ncol(d), "cols\n")

vars <- c("Ta_canopy", "RH_canopy", "WS_canopy", "Rs", "Rln_in", "Prcp")

# 合理值域（参照 Report §4.2 与 BEPS 物理范围）
# Ta:     -90 ~ 60 °C
# RH:      0 ~ 100 %（含饱和 100，雾/雨时常见）
# WS:      0 ~ 40 m/s（哨兵 >40 m/s 已记录于 Plans/bugs.md）
# Rs:      0 ~ 1400 W/m²（地表太阳常数上限 ≈ 1361，地面 ~ 0~1100）
# Rln_in: 50 ~ 500 W/m²（夜间下行长波典型 200~400；负值/小正值为哨兵）
# Prcp:   0 ~ 100 mm/hr（h_max=431 mm/hr 已被识别为固城异常，正常 0~30）

lo <- c(Ta_canopy = -90, RH_canopy = 0, WS_canopy = 0,
        Rs = 0, Rln_in = 50, Prcp = 0)
hi <- c(Ta_canopy = 60,  RH_canopy = 100, WS_canopy = 40,
        Rs = 1400, Rln_in = 500, Prcp = 100)

# 哨兵阈值：明显超出物理范围即记为 sentinel
# 例如 Ta=-103.8, Uz=357 m/s, RH=-26.7, Rln_in=-146
sentinel_lo <- c(Ta_canopy = -90, RH_canopy = -10, WS_canopy = -50,
                 Rs = -10, Rln_in = -200, Prcp = -10)
sentinel_hi <- c(Ta_canopy = 70,  RH_canopy = 110, WS_canopy = 60,
                 Rs = 2000, Rln_in = 800, Prcp = 200)

per_site <- d[, {
  out <- list(n = .N,
              nyear = length(unique(format(time, "%Y"))))
  for (v in vars) {
    x <- suppressWarnings(as.numeric(get(v)))
    finite <- x[is.finite(x)]
    out[[paste0(v, "_min")]]  <- if (length(finite)) min(finite) else NA_real_
    out[[paste0(v, "_max")]]  <- if (length(finite)) max(finite) else NA_real_
    out[[paste0(v, "_mean")]] <- if (length(finite)) mean(finite) else NA_real_
    out[[paste0(v, "_p99")]]  <- if (length(finite)) quantile(finite, 0.99, names = FALSE) else NA_real_
    out[[paste0(v, "_NA_pct")]] <- round(100 * mean(!is.finite(x)), 2)
    # 哨兵：超出物理范围的样本数
    n_sentinel <- sum(!is.na(x) & (x < sentinel_lo[v] | x > sentinel_hi[v]))
    out[[paste0(v, "_sentinel_n")]] <- n_sentinel
    out[[paste0(v, "_sentinel_pct")]] <- round(100 * n_sentinel / .N, 3)
    # 范围外（仍在哨兵内但超合理）：用 hi 阈值
    n_outofrange <- sum(!is.na(x) & (x < lo[v] | x > hi[v]))
    out[[paste0(v, "_outofrange_pct")]] <- round(100 * n_outofrange / .N, 3)
  }
  out
}, by = site]

setorder(per_site, site)
out_file <- "/mnt/z/GitHub/cug-hydro/ModelDev/BEPS.jl/Project_ChinaFlux/OUTPUT/hourly_value_check.csv"
fwrite(per_site, out_file)
cat("Wrote", out_file, "\n")
cat("Rows:", nrow(per_site), "Sites\n\n")

# 摘要：列出每个变量超合理范围的站
cat("=== 摘要：超合理值域 (outofrange > 0.5%) ===\n")
for (v in vars) {
  pct_col <- paste0(v, "_outofrange_pct")
  bad <- per_site[get(pct_col) > 0.5, .(site, get(pct_col))]
  if (nrow(bad) > 0) {
    cat(sprintf("\n[%s] 超合理范围的站点（%d 站）:\n", v, nrow(bad)))
    print(bad)
  }
}

cat("\n=== 摘要：哨兵值 (sentinel > 0.05%) ===\n")
for (v in vars) {
  pct_col <- paste0(v, "_sentinel_pct")
  bad <- per_site[get(pct_col) > 0.05, .(site, get(pct_col))]
  if (nrow(bad) > 0) {
    cat(sprintf("\n[%s] 哨兵值异常的站点（%d 站）:\n", v, nrow(bad)))
    print(bad)
  }
}
