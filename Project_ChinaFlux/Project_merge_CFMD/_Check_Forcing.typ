// _Check_Forcing.typ
// hourly 强迫数据（Forcing_Hourly_Met_sp31_v20260614.csv）的逐站逐变量值域核验 + 异常清洗
// 由 Project_ChinaFlux/check_hourly_range.R 与 sanitize_forcing.R 产出本文件所需全部数据
// 报告生成时间：2026-06-29
#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, size:11pt, footer: "CUG水文气象学2026", header: "")


= 1 数据与方法

*数据源*

`/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv`

1,593,487 行 × 8 列（site, time, Ta_canopy, RH_canopy, WS_canopy, Rs, Rln_in, Prcp）；31 站。
文件版本：`v20260614`（含 2026-06-14 后固城 Prcp ÷29.1、盘锦/长岭 Prcp ÷2.0 的校正）。

#h(2em)

*核验脚本*：`Project_ChinaFlux/check_hourly_range.R`

输出：`Project_ChinaFlux/OUTPUT/hourly_value_check.csv`（31 站 × 47 列，每站每变量含
min/max/mean/p99/NA%/outlier_n/outlier%/outofrange%）。

#h(2em)

*清洗脚本*：`Project_ChinaFlux/sanitize_forcing.R`

运行：`Rscript Project_ChinaFlux/sanitize_forcing.R`

输入：`Forcing_Hourly_Met_sp31_v20260614.csv`

输出：`Project_ChinaFlux/Forcing_Hourly_Met_sp31_v20260614_sanitized.csv`

清洗规则（物理异常 → NA；夜间传感器噪声 → 钳零）：

#table(
  columns: (2.4cm, 2.6cm, 4fr, 3.6cm),
  align: (left, left, left, center),
  table.header([*变量*], [*阈值*], [*策略*], [*清洗量*]),
  [Ta_canopy], [≤ -70 或 ≥ 60 ℃], [异常 → NA；保留 [-70, 60] 物理区间], [649 行 NA],
  [RH_canopy], [< 0 或 > 100 %], [异常 → NA；饱和 100% 物理合理], [466 行 NA],
  [WS_canopy], [< 0 或 > 50 m/s], [异常 → NA；若尔盖/达茂 357/359 m/s 极端异常归此], [103,900 行 NA],
  [Rs], [< 0 或 > 1400 W/m²], [负值钳零（夜间传感器噪声）；> 1400 → NA], [149,003 钳零 + 413 NA],
  [Rln_in], [< 50 或 > 500 W/m²], [异常 → NA；拦截盘锦/长岭/固城/锦州单位错], [70,952 行 NA],
  [Prcp], [不动], [无异常；小时累计 ≤ 100 mm/hr 物理合理], [—],
)

= 2 结论的推导逻辑

第一步：物理合理边界（来自 BEPS 模型输入约束与气象观测常识）：
- Ta ∈ [-70, 60] ℃（地表温度历史极值）
- RH ∈ [0, 100] %（饱和 = 100%）
- WS ∈ [0, 50] m/s（强对流上限；正常通量塔 ≤ 40）
- Rs ∈ [0, 1400] W/m²（地表太阳常数上限 1361）
- Rln_in ∈ [50, 500] W/m²（夜间下行长波典型 200~400）
- Prcp ∈ [0, 250] mm/hr（极端暴雨）

#h(2em)

第二步：异常阈值（明显超出物理范围即判为仪器故障或单位错）：
与上同，但额外标记为 outlier。`outofrange > 0.5%` 或 `outlier > 0.05%` 触发报告。

#h(2em)

第三步：4 类缺陷的成因归类：
- *(a)* 整列 NA = 上游未提供数据（如临泽/三江源/若尔盖 Rs、海北高寒草甸 WS）
- *(b)* 异常值 = 单点仪器故障或单位错混入（如盘锦 WS -2800、固城 Ta -107）
- *(c)* Rs 负值 = 夜间传感器噪声，钳零即可（不视为错误，仅规范）
- *(d)* Prcp_h/Prcp_d 比值异常 = 单位/时间步错（不在本核验范围，见 `TODO.md §A`）

#h(2em)

*注*：Rln_in（下行长波辐射）多数站点未配 pyrgeometer，整列 NA 属正常配置，
不在缺陷清单内。仅当某站 Rln_in 提供了观测值但数值异常（单位错/异常）时才列入 §3.2。

= 3 站点问题清单

== 3.1 全列缺失（NA = 100%，需 CMFD 或他源填补）

#h(2em)
*注*：仅 Rln_in 整列 NA 不列入（多数站未配 pyrgeometer，属正常配置）。下表仅列出
缺 Rs / WS / Prcp 的站（这些是 BEPS 必输入项）。

#table(
  columns: (4.6cm, 3.4cm, 1fr),
  align: (left, left, left),
  table.header([*站点*], [*全 NA 列*], [*影响*]),
  [CRO\_制种玉米\_临泽], [Rs（外加 Prcp=0）], [无辐射 + 无降水，ET_sim=150 vs ET_obs=716 mm/yr，模拟无意义。],
  [GRA\_人工垂穗披碱草\_三江源], [Rs, WS], [无辐射 + 无风，ET_sim=8 vs ET_obs=483 mm/yr。],
  [GRA\_高寒草甸\_海北], [WS], [无风，阻抗计算异常；Rs 存在。],
  [GRA\_高寒草甸\_若尔盖], [Rs], [无辐射，外加 WS 异常 11.39%、Ta 异常 0.79%、RH 异常 0.91%。],
)

== 3.2 单位错（Rln_in 提供观测但值域异常，需溯源数据方）

#h(2em)
本表列出 #emph[有 Rln_in 观测值但数值异常] 的站；Rln_in 整列 NA 的站属正常配置，不列入。

#table(
  columns: (4.6cm, 3cm, 3cm, 1fr),
  align: (left, left, left, left),
  table.header([*站点*], [*Rln_in 值域 (W/m²)*], [*p99*], [*推测单位 / 修复*]),
  [CRO\_水稻\_盘锦], [$-$249.5 ~ 61.4], [24.6], [单位错为最严重之一；正常应 100–500。清洗后 99.96% NA。],
  [CRO\_水稻\_长岭], [$-$48.5 ~ 71.3], [25.0], [同盘锦类问题；正常应 100–500。清洗后 99.91% NA。],
  [CRO\_冬小麦夏玉米\_固城], [$-$146.8 ~ 3.2], [0.0], [符号/单位错；正常应 200–400。清洗后 100% NA。],
  [CRO\_春玉米\_锦州], [$-$683.6 ~ 473.4], [438.4], [单行 $-$999 异常混入；正常段 [100, 473]。清洗后 6.45% NA。],
  [ENF\_人工针叶林\_燕山 ×2], [0 ~ 13.4], [3.5], [98% NA；非 NA 值 < 50，疑 0 异常混入。],
)

== 3.3 异常值（单点仪器故障，清洗后归 NA）

#table(
  columns: (4.6cm, 2.2cm, 1fr, 1fr),
  align: (left, left, left, center),
  table.header([*站点*], [*变量*], [*异常值*], [*样本占比*]),
  [GRA\_高寒草甸\_若尔盖], [WS], [max = 357.6 m/s（与 `Plans/bugs.md #1` 完全吻合）], [10.5%],
  [SAV\_荒漠草原\_达茂], [WS], [max = 359.0 m/s（与若尔盖同量级）], [2.2%],
  [CRO\_水稻\_盘锦], [WS], [min = $-$2800 m/s（疑似 $-99 \u{00D7} 30$ 脏数据）], [3.2%],
  [CRO\_冬小麦夏玉米\_固城], [Ta/RH/WS], [-107.3℃ / -76% / -50 m/s], [0.04% / 0.04% / 0.87%],
  [GRA\_高寒草甸\_若尔盖], [Ta/RH], [-103.8℃ / -26.7%], [0.79% / 0.91%],
  [SAV\_荒漠草原\_达茂], [Ta], [max = 95℃], [0.85%],
  [EBF\_热带雨林\_西双版纳], [Ta/RH], [max = 129℃ / 170%（仪器故障）], [0.01% / 0.01%],
  [CRO\_春玉米\_锦州], [RH], [min = -7.0%], [0.25%],
  [ENF\_北方林森林\_呼中], [RH], [max = 100.6%（湿度计饱和漂移）], [4.94%],
  [CRO\_水稻\_句容], [Rs], [max = 2344 W/m²（远超 1400 上限）], [0.07%],
)

== 3.4 噪声（Rs 负值 = 夜间传感器噪声，已钳零）

- ENF\_北方林森林\_呼中（min = -7.24，46% 负值）
- GRA\_典型草原\_多伦（min = -24.01，45% 负值）
- GRA\_刈割草原\_锡林浩特（min = -5.60，46% 负值）
- GRA\_高寒草甸\_海北（min = -8.12，46% 负值）
- CRO\_冬小麦夏玉米\_固城（min = -5.97，42% 负值）
- CRO\_水稻\_句容（min = -0.99，25% 负值）
- EBF\_亚热带常绿阔叶林\_金佛山（min = -0.10，5.9% 负值）
- CRO\_水稻\_长岭（min = -2.29，0.9% 负值）
- WSA\_高寒灌丛\_海北（min = -0.88，0.6% 负值）

= 4 清洗前后对比

#table(
  columns: (3.6cm, 1fr, 1fr),
  align: (left, center, center),
  table.header([*变量*], [*清洗前 NA%*], [*清洗后 NA%*]),
  [Ta_canopy], [0.00%], [0.06%],
  [RH_canopy], [0.00%], [0.28%],
  [WS_canopy], [0.00%], [6.52%],
  [Rs], [0.00%], [8.07%（其中 9.34% 来自钳零，0.03% 来自 >1400 NA 化）],
  [Rln_in], [61.00%（主因：19 站未配 pyrgeometer，属正常配置）], [72.21%（+11.21% 来自 5 站单位错 NA 化）],
  [Prcp], [0.00%], [0.00%],
)

= 5 优先修复清单

+ *Rln_in 单位错溯源*：盘锦 / 长岭 / 固城 / 锦州 / 燕山×2 共 5 站（注意不是 6），下行长波单位/符号异常。
  已通过 sanitize 拦截为 NA。*注*：其余 19 站 Rln_in 整列 NA 属正常配置（未配 pyrgeometer），不是缺陷。

+ *Rs（短波辐射）全 NA*：临泽 / 三江源 / 若尔盖 共 3 站，
  已加入 SITES_bad 候选，CMFD 实验（`Report_China_FluxALL.typ §6`）显示 Flux_NSE $-$0.6 → $+$0.6 改善。

+ *Prcp 单位 / 时间步错*（不属于本值域核验）：盘锦 ÷2.0、长岭 ÷2.0、固城 ÷29.1 已落实；
  *句容 ÷4.5、金佛山 ÷3.0 尚未落到 v20260614*，见 `Project_ChinaFlux/Plan/TODO.md §A`。

+ *风速异常溯源*：若尔盖 / 达茂 WS 359 m/s 量级相同，疑同一仪器问题；
  盘锦 WS $-$2800 m/s 疑 $-$99×30 脏数据。已 NA 化，溯源后再决定是否补回。

+ *元数据 VegType 误标*：若尔盖 / 三江源 CRO → GRA（高寒草甸），
  见 `Report_China_FluxALL.typ §4.2`。

= 6 输出文件

- 核验结果：`Project_ChinaFlux/OUTPUT/hourly_value_check.csv`
- 清洗后 forcing：`Project_ChinaFlux/Forcing_Hourly_Met_sp31_v20260614_sanitized.csv`
- 清洗后核验：`Project_ChinaFlux/OUTPUT/hourly_value_check_sanitized.csv`
- 核验脚本：`Project_ChinaFlux/check_hourly_range.R`
- 清洗脚本：`Project_ChinaFlux/sanitize_forcing.R`
- 主报告引用：`Project_ChinaFlux/Plan/Report_China_FluxALL.typ §4.2, §4.2.1, §5, §7`
