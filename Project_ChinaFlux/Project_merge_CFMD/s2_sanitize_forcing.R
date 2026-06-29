#!/usr/bin/env Rscript
# sanitize_forcing.R — hourly forcing 值域清洗
# 输入：/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv
# 输出：.../Forcing_Hourly_Met_sp31_v20260614_sanitized.csv
#
# 策略：
#   1. 异常值 → NA（不静默截断，避免下游误判）
#   2. Rs < 0 → 钳为 0（夜间传感器噪声，物理意义为零）
#   3. 保留原有的 NA；不填补
#
# 阈值来自 OUTPUT/hourly_value_check.csv 核对结果（2026-06-29）
pacman::p_load(data.table, Ipaper)

src <- "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv"
dst <- "/mnt/z/GitHub/cug-hydro/ModelDev/BEPS.jl/Project_ChinaFlux/Project_merge_CFMD/Forcing_Hourly_Met_sp31_v20260614_sanitized.csv"

d <- fread(src)

# styler: off
rules <- tribble(
  ~var,        ~lower, ~upper, ~action,        ~note,
  "Ta_canopy",    -70,     60, "to_na",        "高寒站点历史最低 ~-50；地表上限 60",
  "RH_canopy",      0,    100, "to_na",        "RH<0 仪器错；>100 物理不可能",
  "WS_canopy",      0,     50, "to_na",        "<0 测量错（盘锦 -2800）；>50 传感器故障（若尔盖/达茂 357）",
  "Rs",             0,   1400, "clamp_low_na", "<0 钳为 0（夜间噪声）；>1400 NA（地表上限 ~1361）",
  "Rln_in",        50,    500, "to_na",        "<50 单位/符号错；>500 仪器故障",
  "Prcp",        -Inf,    Inf, "passthrough",  "小时累计 ≤100 物理合理"
)
# styler: on

# 复用入口：action 决定三种清洗语义
#   to_na         越界值 → NA（含两端）
#   clamp_low_na  下端钳零、上端 NA
#   passthrough   不动
apply_rule <- function(d, var, lower, upper, action) {
  x <- d[[var]]
  switch(action,
    to_na = {
      n_bad <- sum(!is.finite(x) | x < lower | x > upper)
      d[get(var) < lower | get(var) > upper, (var) := NA_real_]
      cat(sprintf("  [%s] %s: 越界 → NA: %d 行\n", var, action, n_bad))
    },
    clamp_low_na = {
      n_low <- sum(x < lower, na.rm = TRUE)
      n_high <- sum(x > upper, na.rm = TRUE)
      d[get(var) < lower, (var) := 0]
      d[get(var) > upper, (var) := NA_real_]
      cat(sprintf(
        "  [%s] %s: <%g 钳零 %d 行；>%g NA %d 行\n",
        var, action, lower, n_low, upper, n_high
      ))
    },
    passthrough = {
      cat(sprintf("  [%s] passthrough（不动）\n", var))
    }
  )
  invisible(d)
}

for (i in seq_len(nrow(rules))) {
  apply_rule(d, rules$var[i], rules$lower[i], rules$upper[i], rules$action[i])
}
fwrite(d, dst)

# 验证：每列 NA 占比与极值
cat("\n=== 清洗后核对 ===\n")
d2 <- fread(dst)

for (v in rules$var) {
  cat(sprintf(
    "  %s: NA_pct=%.2f, min=%.2f, max=%.2f\n", v,
    100 * mean(is.na(d2[[v]])),
    min(d2[[v]], na.rm = TRUE),
    max(d2[[v]], na.rm = TRUE)
  ))
}
