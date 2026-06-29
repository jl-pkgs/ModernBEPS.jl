# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 包信息

- **注册名**：`BEPS`（`Project.toml` 中的 `name`）；仓库与 CI/文档仍以 `ModernBEPS.jl` 名义发布
- **开发文档**：<https://jl-pkgs.github.io/ModernBEPS.jl/dev>
- **科学定位**：耦合土壤–植物–大气（SPAC）的陆面过程模型，逐小时步长模拟碳、水、热通量
- **双实现路径**：
  - `inter_prg_jl`：纯 Julia 实现（主路径，性能约为 C 版本 4.3×，见 `README.md`）
  - `inter_prg_c` / `beps_main(version="c")`：通过 `libbeps` 共享库调用 C 版本（保留用于一致性校验）

## 项目结构

```
src/        物理过程模块（辐射/光合/SPAC/雪/土壤水热）
  DataType/ 参数与状态类型（ParamBEPS、StateBEPS、Met、Flux、ETFlux、LeafCache）
  SPAC/     冠层–叶片辅助函数（VCmax、LAI 划分、雪密度）
  SoilPhysics/ 土壤水热更新
  clang/    libbeps 的 C 绑定（Clang.jl 生成）
  standalone/Photosynthesis/  独立 Farquhar 模型（无冠层耦合，便于单元测试）
test/       测试；主入口 test/runtests.jl
examples/   入口示例：example_01.jl（最简流程）、Figure1_compare_with_C.jl（C/Julia 对比）
data/       查找表：SOIL_USDA2BEPS.csv、VEG_IGBG2BEPS.csv、Tgound_CUG.csv
docs/       手册（manual.typ）+ 各章节笔记 + 图件
Plans/      调试与验证计划（含 bugs.md 历史记录）
Project_ChinaFlux/  ChinaFlux 通量站验证子项目（独立 case01/02/03 驱动脚本）
```

## 常用命令

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'  # 装依赖
julia --project test/runtests.jl                    # 跑全部测试
julia --project test/test-beps_modern.jl            # 单文件测试（test/ 下任一文件同此格式）
julia --project                                    # 进入项目 REPL
julia --project -e 'using Pkg; Pkg.test()'         # 触发 Pkg 测试（含依赖预编译）
```

> 用户机器上 `julia` 别名指向 `/home/kong/.local/bin/julia`，`Rscript` → `/opt/miniforge3/envs/r4.5/bin/Rscript`。

## 架构与数据流

### 入口与主循环

1. `setup_model(VegType, SoilType)`（或 `setup`）：根据植被/土壤代码从 `data/` 查表，构建 `ParamBEPS{FT,N,H,T}` + `StateBEPS`
2. `simulate(forcing, lai, dates; ps, state, ...)`：**现代规范 API**，逐时步调用 `inter_prg_jl`，返回 `(df_flux, df_ET, states, caches)`
3. `inter_prg_jl`：单小时步长，顺序执行 辐射 → 空气动力学 → 冠层能量平衡迭代 → 光合 → 蒸散 → 雪/截留 → 土壤水热
4. 土壤水热：`UpdateSoilMoisture!` 委托 `SoilDiffEqs.jl`，`UpdateHeatFlux!` 处理热扩散
5. `simulate` 的 `SM_obs` / `TS_obs` kwarg 可旁路土壤水/热更新（用于观测同化或诊断）

### 现代 API（首选）

```julia
state, ps = setup_model("evergreen_needleleaf", "silty_clay_loam")
df_flux, df_ET, states, caches = simulate(forcing, lai, dates;
  ps, state, lon=120.0, lat=40.0,
  SM_obs=nothing, TS_obs=nothing,
  VARS_STATE=[:θ, :Tsoil_c, :z_snow],
  VARS_CACHE=[:Gs_o, :Gs_u])
```

### 遗留 API（仅用于 C/Julia 一致性校验）

`beps_main(forcing, lai, dates; ...)` 在源码中已标注 "Please use `simulate` instead"，仅在 `test-beps_main.jl` 中用于 C 版本对比。新代码不应调用。

### 关键类型

| 类型 | 角色 |
| --- | --- |
| `ParamBEPS{FT,N,H,T}` | 顶层参数容器（`hydraulic::HydraulicProfile` + `thermal::ThermalProfile` + `veg::ParamVeg`） |
| `StateBEPS` | 可变状态：土壤 θ/T、积雪、根系分布、水分胁迫 |
| `ParamVeg{FT}` | 植被参数：LAI、反照率、气孔斜率、`VCmax25`、聚集指数 |
| `MetSeries` / `Met` | 逐小时强迫 vs 单时步结构（`fill_met!` 拷贝） |
| `Flux` / `ETFlux` | 通量输出 |
| `LeafCache` | 每叶中间量缓存（温度、导度）；`CacheSeries` 时间序列 |
| `StateSeries` | 选中状态变量的时间序列容器 |

## 参数管理约定（ModelParams.jl）

**所有参数**用 `@bounds @with_kw` 定义，含默认值与边界 `(min, max)`：

```julia
@bounds @with_kw mutable struct ParamBEPS{FT<:AbstractFloat}
  N::Int = 5
  dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]
  r_drainage::FT = 0.50 | (0.2, 0.7)
  # ...
end
```

- 默认值与边界统一在 `src/DataType/Params/` 内管理
- 物理量单位必须就地标注（`# [m]`, `[W m-1 K-1]` 等），单位不统一是常见 bug 来源（见 `docs/Figures/Figure3_变量单位统一.png`）
- 结构体应力求类型稳定；`FT`/`N`/`H`/`T` 等参数化类型保持值类型参数，避免 `@with_kw` 在值类型上引入 Union boxing

## 测试

主入口 `test/runtests.jl` 按顺序 `include` 11 个文件：

| 文件 | 覆盖范围 |
| --- | --- |
| `test-beps_main.jl` | C/Julia 一致性（遗留 API） |
| `test-macro.jl` | `@bounds` / `@with_kw` 等宏 |
| `test-StateSeries.jl` | 状态序列容器 |
| `dev/test-aerodynamic_conductance.jl` | 空气动力学导度（V1/V2） |
| `test-photosynthesis_standalone.jl` | 独立 Farquhar 模型 |
| `test-beps_modern.jl` | 端到端 `simulate` |
| `test-soilwater.jl` | 土壤水模块 |
| `test-ModelParams.jl` | ModelParams 集成 |
| `test-utilize.jl` | 工具函数 |
| `test-UpdateSoilMoisture_Q0.jl` | 土壤水求解器 |
| `modules/modules.jl` | 模块加载与导出 |

`test/debug/`、`test/dev/`、`test/modules/` 内为诊断/开发用脚本。

## 已知坑与近期修复

`README.md` 末尾「Bugs Fixed」是权威历史记录；`Plans/bugs.md` 收录带推导过程的诊断。最近（v0.1.x）已修：

- `inter_prg` 中 `UpdateHeatFlux(state, Ta_annual, kstep)` 第二参数必须传 `Ta_annual` 而非 `Tair`
- `surface_temperature_jl` `T_weighted` 分子漏 `z_snow`、`G_soil` 漏 `κ_soil1`
- `Init_Soil_Parameters` 中 `V_SOM` 量纲为 [0,1] 而非 0–100
- Obukhov 长度公式分子分母颠倒（`Plans/bugs.md #1`）

未关闭的 TODO 见 `README.md` 顶部「TODO」与 `Project_ChinaFlux/Plan/SITES_bad.md`。

## 风险与注意

- **`docs/manual/manual.typ` 是 Claude 自动生成的草稿**，文档顶部已自标注需谨慎使用；不要将其内容当作权威推导
- `Project.toml` 注册名 `BEPS` 与 GitHub 仓库/文档站 `ModernBEPS.jl` 不一致，导入/引用时分别用
- 修改物理公式前先查 `Plans/bugs.md` + `git log -- src/<file>.jl`，避免重蹈已修问题
- `Project_ChinaFlux/` 是独立验证子项目，调试通量站问题时优先看其 `Plan/` 与 `SITES_bad.md`

## 沟通与代码风格（来自 `.github/instructions/main.instructions.md`）

- **正文使用中文**，代码标识符、变量名、文件路径保留英文
- **Linux 极简主义**：不啰嗦、不堆砌；注释只保留"未来复用"或"帮助理解"两类
- 代码排版规范、可读性优先（用户 20 年码龄，对"行数少"无偏好）
- 物理量单位就地标注；类型稳定优先于微优化