#!/usr/bin/env Rscript
# sanitize_forcing.R — hourly forcing 值域清洗
# 输入：/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv
# 输出：.../Forcing_Hourly_Met_sp31_v20260614_sanitized.csv
#
# 策略：
#   1. 物理哨兵 → NA（不静默截断，避免下游误判）
#   2. 合理范围外的值 → NA
#   3. Rs < 0 → 钳为 0（夜间传感器噪声，物理意义为零）
#   4. 保留原有的 NA；不填补
#
# 阈值来自 OUTPUT/hourly_value_check.csv 核对结果（2026-06-29）

suppressPackageStartupMessages(library(data.table))

src <- "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv"
dst <- "/mnt/z/GitHub/cug-hydro/ModelDev/BEPS.jl/Project_ChinaFlux/Forcing_Hourly_Met_sp31_v20260614_sanitized.csv"

cat("[1/4] 读入", src, "\n")
d <- fread(src)
cat("  rows =", nrow(d), "sites =", uniqueN(d$site), "\n\n")

# --- 1. Ta_canopy 哨兵 + 范围 ---
#   Ta < -70 必是哨兵（高寒站点历史最低 ~-50）
#   Ta > 60 必是哨兵（地表温度上限）
cat("[2/4] 清洗 Ta_canopy\n")
n_before <- sum(!is.finite(d$Ta_canopy) | d$Ta_canopy < -70 | d$Ta_canopy > 60)
d[Ta_canopy < -70 | Ta_canopy > 60, Ta_canopy := NA_real_]
cat("  哨兵/越界 → NA: ", n_before, " 行\n")

# --- 2. RH_canopy 哨兵 + 范围 ---
#   RH < 0 必是哨兵；RH > 100 物理不可能（除非仪器故障）
cat("[3/4] 清洗 RH_canopy\n")
n_before <- sum(!is.finite(d$RH_canopy) | d$RH_canopy < 0 | d$RH_canopy > 100)
d[RH_canopy < 0 | RH_canopy > 100, RH_canopy := NA_real_]
n_after <- sum(!is.finite(d$RH_canopy))
cat("  哨兵 → NA: ", n_before, " 行\n")

# --- 3. WS_canopy 哨兵 + 范围 ---
cat("[4/4] 清洗 WS_canopy\n")
# 哨兵阈值：
#   WS < 0 是测量错误（如盘锦 -2800）
#   WS > 50 m/s 是传感器故障（若尔盖 357、达茂 359）
#   正常范围 0 ~ 40，山区强对流可达 ~50
n_before <- sum(!is.finite(d$WS_canopy) | d$WS_canopy < 0 | d$WS_canopy > 50)
d[WS_canopy < 0 | WS_canopy > 50, WS_canopy := NA_real_]
cat("  哨兵 → NA: ", n_before, " 行\n")

# --- 4. Rs 哨兵 + 范围 ---
#   负值 → 钳为 0（夜间传感器噪声）
#   > 1400 W/m² → NA（地表太阳常数上限 ~1361）
n_neg <- sum(d$Rs < 0, na.rm = TRUE)
n_high <- sum(d$Rs > 1400, na.rm = TRUE)
d[Rs < 0, Rs := 0]
d[Rs > 1400, Rs := NA_real_]
cat("[5/7] 清洗 Rs\n")
cat("  Rs<0 钳零:", n_neg, "行\n")
cat("  Rs>1400 → NA:", n_high, "行\n")

# --- 5. Rln_in 哨兵 + 范围（最重要的修复）---
#   Rln_in 是下行长波辐射，正常 100~500 W/m²
#   负值或 < 50 → 单位错或符号错（盘锦、长岭、固城、锦州）
#   > 500 → 仪器故障
n_neg_rln <- sum(d$Rln_in < 50, na.rm = TRUE)
n_high_rln <- sum(d$Rln_in > 500, na.rm = TRUE)
d[Rln_in < 50, Rln_in := NA_real_]
d[Rln_in > 500, Rln_in := NA_real_]
cat("[6/7] 清洗 Rln_in（单位错哨兵）\n")
cat("  Rln_in<50 → NA:", n_neg_rln, "行\n")
cat("  Rln_in>500 → NA:", n_high_rln, "行\n")

# --- 6. Prcp 不动 ---
# Prcp 没有哨兵；小时累计 ≤ 100 mm/hr 合理（极端暴雨可达）

# --- 7. 输出 ---
cat("[7/7] 写出", dst, "\n")
fwrite(d, dst)
cat("完成。\n")

# --- 验证 ---
cat("\n=== 清洗后核对 ===\n")
d2 <- fread(dst)
cat("总行数:", nrow(d2), "\n")
for (v in c("Ta_canopy","RH_canopy","WS_canopy","Rs","Rln_in","Prcp")) {
  cat(sprintf("  %s: NA_pct=%.2f, min=%.2f, max=%.2f\n", v,
    100*mean(is.na(d2[[v]])), min(d2[[v]], na.rm=TRUE), max(d2[[v]], na.rm=TRUE)))
}