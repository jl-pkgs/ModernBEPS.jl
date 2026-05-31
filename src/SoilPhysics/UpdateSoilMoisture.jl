# soil moisture transport is governed by infiltration, runoff, gradient diffusion, 
# gravity, and root
# extraction through canopy transpiration.  The net water applied to the
# surface layer is the snowmelt plus precipitation plus the throughfall
# of canopy dew minus surface runoff and evaporation.
# CLM3.5 uses a zero-flow bottom boundary condition.
include("Modules.jl")
include("SolveSM_BEPS.jl")
include("SolveSM_Bonan.jl")


# !deprecated：兼容 Soil 结构体
UpdateSoilMoisture(soil::Soil, kstep::Float64) = UpdateSoilMoisture(soil, soil, kstep)


"""
    UpdateSoilMoisture(st, ps, kstep; fix_sm=false, SolveSM_fn=SolveSM_BEPS)

更新土壤水分状态 `st.θ`，JAX 风格 `(st, ps)` 签名。

包含两个阶段：
1. **地表积水与入渗**：更新冻结因子、计算入渗率 `inf` [cm/h]、地表径流后的积水 `st.z_water`。
2. **土壤水求解**：调用 `SolveSM_fn(st, ps, inf, kstep)` 求解各层 θ，再更新冰比例。

# 关键字
- `fix_sm=true`：仅更新地表积水，不改变土壤水分状态（观测土壤水模式）。
- `SolveSM_fn`：土壤水求解器，可选 `SolveSM_BEPS`（显式步长自适应）或
  `SolveSM_Bonan`（Bonan-Q0 隐式 Crank-Nicolson）。
"""
function UpdateSoilMoisture(st::S, ps::P, kstep::Float64;
  fix_sm::Bool=false, SolveSM_fn=SolveSM_BEPS) where {
  S<:Union{StateBEPS,Soil},P<:Union{ParamBEPS,Soil}}

  n = st.n_layer
  (; θ_sat, K_sat, b, ψ_sat) = get_hydraulic(ps)
  (; θ, θ_prev, ice_ratio, f_water, Tsoil_c, dz, z_water, r_rain_g) = st
  r_drainage = ps.r_drainage

  θ_prev .= θ

  # ===== 1. 地表积水与入渗 =====
  update_SoilWaterFrac!(f_water, Tsoil_c; n)
  inf = cal_infiltration(θ, dz, K_sat, θ_sat, ψ_sat, b, f_water[1], z_water, r_rain_g, kstep)
  st.inf = inf
  inf_ms = inf / 360000.0 # [cm h-1] -> [m s-1]
  st.z_water = (z_water / kstep + r_rain_g - inf_ms) * kstep * r_drainage # Ponded water after runoff

  fix_sm && return # 如果 fix_sm=true，则只更新地表积水，不改变土壤水分状态

  # ===== 2. 土壤水求解 =====
  SolveSM_fn(st, ps, inf, kstep)

  # update ice ratio
  for i in 1:n
    ice_ratio[i] *= θ_prev[i] / θ[i]
    ice_ratio[i] = min(1.0, ice_ratio[i])
  end
end
