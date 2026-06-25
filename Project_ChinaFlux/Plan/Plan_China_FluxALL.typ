#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2026", header: "")

= 研究目标

在 31 个具备小时尺度气象驱动的涡度相关通量站点上运行并率定 BEPS 模型，
覆盖农田（CRO）、落叶/常绿阔叶林（DBF/EBF）、针叶林（ENF）、混交林（MF）、
草地（GRA）、稀树草原（SAV）、湿地（WET）与灌丛（WSA）共九类植被，
以系统评估模型在中国多种生态系统下对碳、水、能量通量的模拟能力。

= 数据来源

驱动与验证数据说明见
`/mnt/z/GitHub/jl-pkgs/ChinaFlux2026/data/BEPS/BEPS_Forcing_China_FluxALL.md`。

#table(
  columns: (auto, 1fr),
  align: (left, left),
  table.header([*用途*], [*文件*]),
  [小时气象驱动], [`data/BEPS/Forcing_Hourly_Met_sp31_v20260614.csv`（含 Ta、RH、WS、Rs、Rln_in、Prcp）],
  [日尺度通量/LAI 验证], [`data/BEPS/Daily_FluxALL/<站点>_Daily_FluxALL_v20260615.csv`（含 GPP、ET、LAI 及分层 SM/TS）],
  [站点元数据], [`data/Metadata/ChinaFlux_Metadata.csv`（经纬度、植被与土壤类型、观测高度、土壤水分/温度深度）],
)

时间约定：模型内部采用 *UTC*，精度评价与输出对齐采用 *当地时间*（UTC+8），
即 `dates_UTC = dates_local - 8h`。

= 执行流程

批量处理脚本为 `Project_ChinaFlux/case02_ChinaFlux_ALL.jl`，由森林专用版
`case01_ChinaFlux_Forest.jl` 扩展而来。单站处理链路为：
*读取驱动* → *初始化参数与状态* → *逐小时模拟* → *日尺度精度评价* → *参数优化*。

针对全站点数据的异构性，脚本在数据接入环节做了鲁棒化处理：

- 缺失通量列（如部分站点无 GPP）自动以 `NaN` 占位，仅评价可用变量；
- 土壤水分/温度按「标准列存在且含有效观测」的深度自动筛选，无匹配则仅评价 GPP/ET/Hs；
- 元数据中字符串型或 `NA` 数值统一转为数值，缺失项以合理默认值兜底；
- 驱动与观测起止时间不一致时，按观测首日对齐起点并裁剪长度。

批量运行采用逐站 `try/catch`：单站报错时记录站点与错误信息后继续，
全部完成后汇总为 `_errors.csv`，便于定位与修复。

= 当前进展

经数据接入修复与若干代码/数据问题订正后，*31 个站点的逐小时模拟（simulate + GOF）已全部跑通*，
结果存于 `Project_ChinaFlux/OUTPUT/ALL/sim_only/BEPS_<站点>.jld2`。
失败站点的错误根因与修复记录见 `Project_ChinaFlux/Plan/SITES_bad.md`。

= 参数优化

在模拟跑通的基础上，对全部站点执行参数优化（SCE-UA 算法，`maxn = 1000`），
以最大化 GPP 与 ET 的 NSE 均值（`goal = :NSE`，`goal_multiplier = -1`）为目标。
优化变量为：

```julia
[:r_drainage, [:veg, :Ω], [:veg, :g1_w], [:veg, :g0_w], [:veg, :VCmax25]]
```

每站输出最优参数 `theta_opt` 及优化前后的精度（`gof`、`gof_opt`）、
日尺度模拟与观测对齐结果（`data_sim`、`data_obs`），
保存至 `Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE/BEPS_<站点>.jld2`。

= 并行运行

本地计算资源充足，可并行提交任务：

```bash
mpiexecjl -n 6 julia -t5 --project Project_ChinaFlux/case02_ChinaFlux_ALL.jl
```
