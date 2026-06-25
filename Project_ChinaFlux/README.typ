#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#align(center)[
  #text(17pt, weight: "bold")[BEPS 模型输入与输出说明]

  #text(10pt)[基于 `Project_ChinaFlux/case01_ChinaFlux.jl` 的中国通量站点率定流程]
]

= 概述

`case01_ChinaFlux.jl` 以中国通量观测站点为例，演示 BEPS 模型的完整运行链路：
*读取驱动* → *初始化参数与状态* → *逐小时模拟* → *日尺度精度评价* → *参数优化*。
核心调用为：

```julia
model = ParamBEPS(VegType, SoilType)        # 参数
state = InitState0(model, forcing)          # 初始状态
df_fluxes, df_ET, states, caches =
  simulate(forcing, lai, dates_UTC; ps=model, state, lon, lat, SolveSM_fn)
gof, data_sim, data_obs =
  BEPS_GOF(df_fluxes, states, dates_local, FluxALL; depths_SM, depths_TS)
theta_opt = optim(model, forcing, lai, dates_UTC; paths, maxn, kw_loss...)
```

模型时间步长为 *1 小时*。注意：`simulate` 内部使用 *UTC* 时间（`dates_UTC = dates_local - 8h`），
而精度评价与输出对齐使用 *当地时间*。

#line(length: 100%, stroke: 0.5pt + gray)

= 模型输入

== 气象驱动 `MetSeries`

逐小时气象强迫，由 `MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)` 构造。
每个字段均为长度 `ntime` 的向量，逐时刻通过 `fill_met!` 写入单步 `Met`。

#table(
  columns: (auto, auto, auto, 1fr),
  align: (left, center, center, left),
  table.header([*字段*], [*单位*], [*必需*], [*说明*]),
  [`Rs`],     [W m⁻²],  [是], [向下短波辐射],
  [`Tair`],   [°C],     [是], [2 m 气温],
  [`RH`],     [%],      [是], [相对湿度（0–100）],
  [`Uz`],     [m s⁻¹],  [是], [测风高度 z 处风速],
  [`Prcp`],   [mm h⁻¹], [是], [降水强度],
  [`Rln_in`], [W m⁻²],  [否], [向下长波辐射；缺测（`NaN`）时由模型内部估算],
)

驱动数据在读入时经 `sanitize_forcing!` 质控：`Tair` ∈ [-50, 50]、`RH` ∈ [0, 100]、
`Uz ≥ 0.01`、`Rs ≥ 0`、`Prcp ≥ 0`（无效值插值或置零）。

== 叶面积指数 `lai`

*日尺度* LAI 向量，长度等于天数；`simulate` 内部按 `k = ceil(i/24)` 将日值广播到逐小时。
案例中取自 GLASS 产品列 `LAI_glass_G005`。

== 时间 `dates`

`DateTime` 向量，长度与 `forcing` 一致。传入 `simulate` 的须为 *UTC* 时间。

== 站点信息（`st_full` 表）

按站点名检索的元数据，用于配置模型几何与土层：

#table(
  columns: (auto, 1fr),
  align: (left, left),
  table.header([*字段*], [*说明*]),
  [`lon`, `lat`],   [经纬度 [°]，用于太阳几何计算],
  [`VegType`],      [植被类型，决定 `ParamVeg` 默认值],
  [`SoilType`],     [土壤质地，决定水力 / 热力参数],
  [`z_Uz`],         [风速测量高度 → `model.veg.z_wind`],
  [`z_overstory`],  [上层冠层高度 → `model.veg.z_canopy_o`],
  [`z_SM`],         [土壤水分观测深度（cm，逗号分隔）→ `depths_SM`],
  [`z_TS`],         [土壤温度观测深度（cm，逗号分隔）→ `depths_TS`],
)

== 模型参数 `ParamBEPS`

由 `ParamBEPS(VegType, SoilType)` 从查找表构造，包含土壤水力 / 热力廓线与植被参数。
顶层及可率定的关键参数：

#table(
  columns: (auto, auto, auto, 1fr),
  align: (left, center, center, left),
  table.header([*参数*], [*默认/范围*], [*单位*], [*含义*]),
  [`r_drainage`],     [0.5 / (0.2, 0.7)],   [—],          [地表排水（汇流）比例],
  [`ψ_min`],          [3300 (1000 阔叶林)], [cm],         [水分胁迫起始点],
  [`alpha`],          [0.4 (1.5 阔叶林)],   [—],          [土壤水限制因子参数 (He 2017)],
  [`veg.Ω`],          [0.85 / (0.3, 1.0)],  [—],          [冠层聚集度指数 clumping],
  [`veg.g1_w`],       [8 / (1, 20)],        [—],          [Ball-Berry 斜率],
  [`veg.g0_w`],       [0.0175 / (.001,.1)], [—],          [Ball-Berry 截距 (H₂O)],
  [`veg.VCmax25`],    [89.45 / (5, 200)],   [μmol m⁻² s⁻¹], [25 ℃ 最大 Rubisco 羧化速率],
  [`veg.r_root_decay`], [0.95 / (.85,.999)], [—],          [根系分布衰减率],
  [`veg.LAI_max_o/u`], [4.5 / 2.4],         [—],          [上/下层最大 LAI],
  [`veg.α_canopy_*`], [—],                  [—],          [冠层反照率 (vis/nir)],
  [`veg.z_canopy_o/u`], [—],                [m],          [上/下层冠层高度],
)

土壤廓线（`hydraulic`：Campbell 持水模型；`thermal`：热力参数）按 5 层离散，
默认层厚 `dz = (0.05, 0.10, 0.20, 0.40, 1.25)` m。

== 初始状态 `StateBEPS`

由 `InitState0(model, forcing)` 生成：初始土温 `Tsoil0 = Tair[1]`，
初始土壤含水量 `θ0 = θ_sat × 0.8`（田间持水量量级），积雪深度 `z_snow0 = 0`。
`StateBEPS` 为可变状态容器，记录逐层 θ、Tsoil、积雪、根系比例、水分胁迫因子等。

== 率定观测数据 `FluxALL`（日尺度）

参数优化与精度评价所需的观测，逐日：

#table(
  columns: (auto, auto, 1fr),
  align: (left, center, left),
  table.header([*列*], [*单位*], [*说明*]),
  [`lai`],     [—],       [叶面积指数（亦作驱动）],
  [`GPP_obs`], [gC m⁻² d⁻¹], [总初级生产力（负值自动取反归一化）],
  [`ET_obs`],  [mm d⁻¹],  [蒸散发],
  [`Hs_obs`],  [W m⁻²],   [感热通量],
  [`SM_*cm`],  [m³ m⁻³],  [各深度土壤体积含水量],
  [`TS_*cm`],  [°C],      [各深度土壤温度],
)

#line(length: 100%, stroke: 0.5pt + gray)

= 模型输出

`simulate` 返回四元组 `(df_fluxes, df_ET, states, caches)`，均为 *逐小时* 结果。

== `df_fluxes`（`Flux` → DataFrame）

碳 / 能量 / 水量通量。

#table(
  columns: (auto, auto, 1fr),
  align: (left, center, left),
  table.header([*列*], [*单位*], [*说明*]),
  [`GPP`],          [gC m⁻²],     [小时总 GPP（= Σgpp × 12 × step × 1e-6）],
  [`NPP`, `NEP`],   [gC m⁻²],     [净初级生产力 / 净生态系统生产力],
  [`gpp_{o,u}_{sunlit,shaded}`], [μmol m⁻² s⁻¹], [上/下层、阳/阴叶瞬时 GPP 分量],
  [`npp_o`, `npp_u`], [gC m⁻²],   [上/下层 NPP],
  [`plant_resp`, `soil_resp`], [gC m⁻²], [植物 / 土壤呼吸],
  [`Net_Rad`],      [W m⁻²],      [净辐射],
  [`SH`],           [W m⁻²],      [感热通量],
  [`LH`],           [W m⁻²],      [潜热通量],
  [`Trans`],        [mm],         [小时蒸腾总量],
  [`Evap`],         [mm],         [小时蒸发总量],
  [`z_water`],      [m],          [地表积水深度],
  [`z_snow`],       [m],          [积雪深度],
  [`ρ_snow`],       [kg m⁻³],     [雪密度],
  [`inf`],          [cm h⁻¹],     [下渗率],
)

== `df_ET`（`ETFlux` → DataFrame）

蒸散发的精细拆分（分量为瞬时速率 `kg m⁻² s⁻¹`，汇总量为小时总量 mm）。

#table(
  columns: (auto, auto, 1fr),
  align: (left, center, left),
  table.header([*列*], [*单位*], [*说明*]),
  [`Trans_o`, `Trans_u`],   [kg m⁻² s⁻¹], [上/下层蒸腾],
  [`Eil_o`, `Eil_u`],       [kg m⁻² s⁻¹], [上/下层液态截留蒸发],
  [`EiS_o`, `EiS_u`],       [kg m⁻² s⁻¹], [上/下层固态（升华）],
  [`Evap_soil`],            [kg m⁻² s⁻¹], [土壤蒸发],
  [`Evap_SW`, `Evap_SS`],   [kg m⁻² s⁻¹], [积水 / 积雪蒸发],
  [`Qhc_o`, `Qhc_u`, `Qhg`],[W m⁻²],      [上层/下层冠层、地面感热],
  [`Trans`, `Evap`],        [mm],         [小时蒸腾 / 蒸发总量],
  [`SH`, `LH`],             [W m⁻²],      [总感热 / 总潜热],
)

== `states`（`StateSeries`）

按 `VARS_STATE`（默认 `[:θ, :Tsoil_c, :z_snow]`）保存的状态时间序列，分两组：
- `states.scalars`：标量字段，`Vector{Float64}`（长度 `ntime`），如 `z_snow`、`z_water`；
- `states.vectors`：逐层字段，`Matrix{Float64}`（`ntime × nlayer`），如 `θ`（土壤含水量）、`Tsoil_c`（土温）。

== `caches`（`CacheSeries`）

按 `VARS_CACHE`（默认 `[:Gs_o, :Gs_u]`）保存的 `LeafCache` 中间量时间序列
（如气孔导度、叶温等），`ntime × 4`（阳/阴 × 上/下叶）。

#line(length: 100%, stroke: 0.5pt + gray)

= 精度评价输出 `BEPS_GOF`

`BEPS_GOF` 将逐小时模拟 `agg_daily` 聚合到 *日尺度*（GPP、ET 求和，Hs 求平均），
按观测深度插值 SM / TS，再与观测对齐计算优度，返回 `(gof, data_sim, data_obs)`。

== `gof`（NamedTuple）

每个分量为含 `KGE`、`NSE`、`R²`、`RMSE` 等指标的 DataFrame：

#table(
  columns: (auto, 1fr),
  align: (left, left),
  table.header([*字段*], [*评价对象*]),
  [`gof.Flux`], [GPP、ET、Hs（优化目标取 `Flux[1:2]` 即 GPP 与 ET 的均值）],
  [`gof.SM`],   [各深度土壤含水量 `SM_*cm`],
  [`gof.TS`],   [各深度土壤温度 `TS_*cm`],
)

== `data_sim` / `data_obs`（DataFrame）

日尺度模拟与观测的对齐结果，列结构一致：
`date | GPP | ET | Hs | SM_*cm... | TS_*cm...`，便于直接绘图比对。

#line(length: 100%, stroke: 0.5pt + gray)

= 参数优化

`optim` 以 SCE-UA 算法最小化损失（默认目标 `goal=:NSE`，`goal_multiplier=-1`，
即最大化 GPP 与 ET 的 NSE 均值），优化变量由 `paths` 指定，案例中为：

```julia
[:r_drainage, [:veg, :Ω], [:veg, :g1_w], [:veg, :g0_w], [:veg, :VCmax25]]
```

返回最优参数向量 `theta_opt`；结果连同 `gof`、`gof_opt`、`data_sim`、`data_obs`
保存为 `OUTPUT/.../BEPS_<SITE>.jld2`。
