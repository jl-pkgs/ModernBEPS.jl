export SnowLand

@with_kw mutable struct SnowLand{FT<:AbstractFloat} <: FieldVector{5,FT}
  T_surf::FT = 0.0      # 雪表面温度, 裸土地表温度
  T_snow0::FT = 0.0     # 雪表温度, !注意是雪表, 不是地表
  T_snow1::FT = 0.0     # 雪层1温度
  T_snow2::FT = 0.0     # 雪层2温度
  T_mix0::FT = 0.0      # !mixed
  # T_soil0::FT = 0.0     # !裸土部分
end

function clamp!(des::SnowLand{FT}, src::SnowLand{FT}, Tair::FT) where {FT<:AbstractFloat}
  lower = Tair - FT(2.0)
  upper = Tair + FT(2.0)

  des.T_surf = clamp(src.T_surf, lower, upper)
  des.T_snow0 = clamp(src.T_snow0, lower, upper)
  des.T_snow1 = clamp(src.T_snow1, lower, upper)
  des.T_snow2 = clamp(src.T_snow2, lower, upper)
  des.T_mix0 = clamp(src.T_mix0, lower, upper)
end
clamp!(land::SnowLand{FT}, Tair::FT) where {FT<:AbstractFloat} = clamp!(land, land, Tair)


import ModelParams: AbstractSoil

# ?     : 需要优化的参数
# state : 状态变量
# //    : 未使用的参数
@with_kw mutable struct Soil <: AbstractSoil{Float64,5}
  flag        ::Cint    = Cint(0) # // not used
  n_layer     ::Cint    = Cint(5) # 土壤层数
  step_period ::Cint    = Cint(1) # // not used

  z_water ::Cdouble = Cdouble(0) # [state]
  z_snow  ::Cdouble = Cdouble(0) # [state]
  inf::Float64 = 0.0             # [cm h-1], 本步下渗率，UpdateSoilMoisture 写入

  # the rainfall rate, un--on understory on ground surface  m/s
  r_rain_g    ::Cdouble = Cdouble(0)        # [state], 达到地地表降水, PE, [m/s]

  soil_r      ::Cdouble = Cdouble(0)        # // not used, soil surface resistance for water
  r_drainage  ::Cdouble = Cdouble(0)        # ? 地表排水速率（地表汇流）
  r_root_decay::Cdouble = Cdouble(0)        # ? 根系分布衰减率, decay_rate_of_root_distribution
  ψ_min       ::Cdouble = Cdouble(0)        # ? 开始胁迫，33[m] = 0.33[MPa]；现按 [cm] 存储
  alpha       ::Cdouble = Cdouble(0)        # ? 土壤水限制因子参数，He 2017 JGR-B, Eq. 4
  f_soilwater ::Cdouble = Cdouble(0)        # [state], 总体的土壤水限制因子

  dz          ::Vector{Float64} = zeros(10) # 土壤厚度
  f_root      ::Vector{Float64} = zeros(10) # [state], 根系比例，root fraction
  w_norm      ::Vector{Float64} = zeros(10) # [state], 每层的土壤水限制因子，已归一化

  κ_dry       ::Vector{Float64} = zeros(10) # ? thermal conductivity
  θ_vfc       ::Vector{Float64} = zeros(10) # ? volumetric field capacity
  θ_res       ::Vector{Float64} = zeros(10) # ? volumetric wilting point
  θ_sat       ::Vector{Float64} = zeros(10) # ? volumetric saturation
  K_sat       ::Vector{Float64} = zeros(10) # ? saturated hydraulic conductivity
  ψ_sat       ::Vector{Float64} = zeros(10) # ? soil matric potential at saturation，现按 [negative cm] 存储
  b           ::Vector{Float64} = zeros(10) # ? Cambell parameter b
  ρ_soil      ::Vector{Float64} = zeros(10) # ? 土壤容重，soil density, for volume heat capacity
  V_SOM       ::Vector{Float64} = zeros(10) # ? 有机质含量，organic matter, for volume heat capacity

  ice_ratio   ::Vector{Float64} = zeros(10) # [state]，ice ratio，
  θ           ::Vector{Float64} = zeros(10) # [state], soil moisture
  θ_prev      ::Vector{Float64} = zeros(10) # [state], soil moisture in previous time
  Tsoil_p     ::Vector{Float64} = zeros(10) # [state], soil temperature in previous time
  Tsoil_c     ::Vector{Float64} = zeros(10) # [state], soil temperature in current time

  f_water     ::Vector{Float64} = zeros(10) # [state], 冻结因子，用于 UpdateSoilMoisture
  ψ           ::Vector{Float64} = zeros(10) # [state], soil matric potential, [negative cm]
  θb          ::Vector{Float64} = zeros(10) # // not used, θ at the bottom of each layer
  ψb          ::Vector{Float64} = zeros(10) # // not used
  r_waterflow ::Vector{Float64} = zeros(10) # [state], vertical water flow rate，现按 [cm h-1] 存储
  Kmid        ::Vector{Float64} = zeros(10) # [state], hydraulic conductivity at middle point, [cm h-1]
  Kb          ::Vector{Float64} = zeros(10) # // not used
  Kavg        ::Vector{Float64} = zeros(10) # [state], average conductivity of two soil layers, [cm h-1]
  Cv          ::Vector{Float64} = zeros(10) # [state], volume heat capacity
  κ           ::Vector{Float64} = zeros(10) # [state]
  ETi         ::Vector{Float64} = zeros(10) # [state], 每层蒸发量ET in each layer
  G           ::Vector{Float64} = zeros(10) # [state], 土壤热通量

  ## temporary variables in soil_water_factor_v2
  f_temp      ::Vector{Float64} = zeros(10) # [state], f_i(Tsoil_i), 温度对水分限制影响, Eq. 5
  w_root      ::Vector{Float64} = zeros(10) # [state], 叠加根系分布比例，f_root[i] * f_stress[i]
  f_stress    ::Vector{Float64} = zeros(10) # [state], f_{w,i}, He et al., 2017, Eq. 3, (水分 + 温度)
end


## 设计哲学: 这里把状态变量与模型参数分隔开
# state, params = setup(model)
# st = StateBEPS, ps = ParamBEPS

# (; N, ibeg,
#   θ, ψ, θ_prev, ψ_prev, ψ_next,
#   ∂θ∂ψ, K, K₊ₕ,
#   a, b, c, d, e, f) = soil

# 拖着`ρ_snow`，`ρ_snow`也是一个状态连续的变量
# https://www.eoas.ubc.ca/courses/atsc113/snow/met_concepts/07-met_concepts/07b-newly-fallen-snow-density/
@with_kw mutable struct StateBEPS <: AbstractSoil{Float64,5}
  n_layer    ::Cint = Cint(5) # 土壤层数
  N          ::Int  = 5       # Bonan-Q0 算法所需，与 n_layer 保持同值
  dz         ::Vector{Float64} = zeros(10) # 土壤厚度（从 ps 复制，方便计算）

  Tsnow_c::SnowLand{FT} = SnowLand{FT}() # [inter_prg], 4:8
  Tsnow_p::SnowLand{FT} = SnowLand{FT}() # [inter_prg], 10:15

  Qhc_o  ::FT = 0.0                      # [inter_prg], [11] sensible heat flux

  m_water::Layer2 = Layer2{FT}()         # [inter_prg], [15, 18] + 1
  m_snow ::Layer3 = Layer3{FT}()         # [inter_prg], [16, 19, 20] + 1
  ρ_snow ::FT = 250.0                    # [inter_prg], [kg m-3] snow density

  z_water    ::Cdouble = Cdouble(0)        # [state]
  z_snow     ::Cdouble = Cdouble(0)        # [state]
  inf        ::Float64 = 0.0               # [cm h-1], 本步下渗率，UpdateSoilMoisture 写入

  # the rainfall rate, un--on understory on ground surface  m/s
  r_rain_g   ::Cdouble = Cdouble(0)        # [state], 达到地地表降水, PE, [m/s]
  f_soilwater::Cdouble = Cdouble(0)        # [state], 总体的土壤水限制因子

  f_root     ::Vector{Float64} = zeros(10) # [state], 根系比例，root fraction

  ice_ratio  ::Vector{Float64} = zeros(10) # [state]，ice ratio，
  θ          ::Vector{Float64} = zeros(10) # [state], soil moisture
  θ_prev     ::Vector{Float64} = zeros(10) # [state], soil moisture in previous time
  Tsoil_p    ::Vector{Float64} = zeros(10) # [state], soil temperature in previous time
  Tsoil_c    ::Vector{Float64} = zeros(10) # [state], soil temperature in current time

  f_water    ::Vector{Float64} = zeros(10) # [state], 冻结因子，用于 UpdateSoilMoisture
  ψ          ::Vector{Float64} = zeros(10) # [state], soil matric potential，现按 [negative cm] 存储
  r_waterflow::Vector{Float64} = zeros(10) # [state], vertical water flow rate，现按 [cm h-1] 存储
  Kmid       ::Vector{Float64} = zeros(10) # [state], hydraulic conductivity at middle point (旧求解器)，现按 [cm h-1] 存储
  Kavg       ::Vector{Float64} = zeros(10) # [state], average conductivity of two soil layers (旧求解器)，现按 [cm h-1] 存储
  Cv         ::Vector{Float64} = zeros(10) # [state], volume heat capacity
  κ          ::Vector{Float64} = zeros(10) # [state]
  ETi        ::Vector{Float64} = zeros(10) # [state], 每层蒸发量ET in each layer
  G          ::Vector{Float64} = zeros(10) # [state], 土壤热通量

  ## temporary variables in soil_water_factor_v2
  f_temp     ::Vector{Float64} = zeros(10) # [state], f_i(Tsoil_i), 温度对水分限制影响, Eq. 5
  w_root     ::Vector{Float64} = zeros(10) # [state], 叠加根系分布比例，f_root[i] * f_stress[i]
  w_norm     ::Vector{Float64} = zeros(10) # [state], 每层的土壤水限制因子，已归一化
  f_stress   ::Vector{Float64} = zeros(10) # [state], f_{w,i}, He et al., 2017, Eq. 3

  # ─── Bonan-Q0 求解器字段 ───────────────────────────────────────────────────
  ibeg     ::Int    = 1         # BEPS 固定从第1层开始
  dt       ::Float64 = 3600.0  # [s] 时间步长（小时）

  Δz_cm    ::Vector{Float64} = zeros(10)  # [cm] 各层厚度，setup 时从 dz*100 填充
  Δz₊ₕ_cm  ::Vector{Float64} = zeros(10)  # [cm] 层心间距，setup 时计算

  ψ_prev   ::Vector{Float64} = zeros(10)  # [cm] 上一半步 ψ
  ψ_next   ::Vector{Float64} = zeros(10)  # [cm] 下一半步 ψ
  ∂θ∂ψ    ::Vector{Float64} = zeros(10)  # [cm⁻¹] specific moisture capacity
  K        ::Vector{Float64} = zeros(10)  # [cm h⁻¹] 层中心水力传导率
  K₊ₕ     ::Vector{Float64} = zeros(10)  # [cm h⁻¹] 层界面水力传导率

  # 三对角矩阵求解临时量
  tri::TriSolver{Float64} = TriSolver{Float64,10}()
end


const VARS_SCALAR = Tuple(
    f for (f, T) in zip(fieldnames(StateBEPS), fieldtypes(StateBEPS))
    if T <: AbstractFloat && f ∉ (:Qhc_o,)
)

const VARS_VECTOR = Tuple(
    f for (f, T) in zip(fieldnames(StateBEPS), fieldtypes(StateBEPS))
    if T == Vector{Float64} && f ∉ (:dz, :f_root, :Δz_cm, :Δz₊ₕ_cm,
                                    :ψ_prev, :ψ_next, :∂θ∂ψ, :a, :b, :c, :d, :e, :f)
)
# :θ_prev, :Tsoil_p
const ALL_VARS_STATE = (VARS_SCALAR..., VARS_VECTOR...)

# 从 Soil 构造 SoilState（兼容旧代码）
function StateBEPS(soil::Soil)
  @unpack n_layer, dz, z_water, z_snow, r_rain_g, f_soilwater,
          f_root, w_norm, ice_ratio, θ, θ_prev, Tsoil_p, Tsoil_c,
          f_water, ψ, r_waterflow, Kmid, Kavg, Cv, κ, ETi, G,
          f_temp, w_root, f_stress = soil

  StateBEPS(;
    n_layer, N=Int(n_layer), dz, z_water, z_snow,
    r_rain_g, f_soilwater,
    f_root, w_norm, ice_ratio, θ, θ_prev,
    Tsoil_p, Tsoil_c, f_water, ψ,
    r_waterflow, Kmid, Kavg, Cv, κ,
    ETi, G, f_temp, w_root, f_stress
  )
end

# 将 StateBEPS 同步回 Soil（兼容旧代码）
function State2Soil!(soil::Soil, st::StateBEPS)
  @unpack z_water, z_snow, r_rain_g, f_soilwater,
          f_root, w_norm, ice_ratio, θ, θ_prev, Tsoil_p, Tsoil_c,
          f_water, ψ, r_waterflow, Kmid, Kavg, Cv, κ, ETi, G,
          f_temp, w_root, f_stress = st

  @pack! soil = z_water, z_snow, r_rain_g, f_soilwater,
                f_root, w_norm, ice_ratio, θ, θ_prev, Tsoil_p, Tsoil_c,
               f_water, ψ, r_waterflow, Kmid, Kavg, Cv, κ, ETi, G,
                f_temp, w_root, f_stress
  return soil
end


# state_hydraulic: 统一 AbstractSoil dispatch 的锚点函数
# StateBEPS 的 Δz 字段名为 dz，通过 NamedTuple 映射为算法所需的 Δz（单位不影响比值计算）
@inline function state_hydraulic(st::StateBEPS)
  (;
    N    = st.N,
    ibeg = st.ibeg,
    dt   = st.dt,
    Δz   = st.dz,          # cal_K! 需要 Δz，dz[m] 与 Δz_cm[cm] 比值等价
    Δz_cm   = st.Δz_cm,
    Δz₊ₕ_cm = st.Δz₊ₕ_cm,
    θ = st.θ, ψ = st.ψ,
    θ_prev = st.θ_prev,
    ψ_prev = st.ψ_prev,
    ψ_next = st.ψ_next,
    K   = st.K,
    K₊ₕ = st.K₊ₕ,
    ∂θ∂ψ = st.∂θ∂ψ,
    a = st.a, b = st.b, c = st.c,
    d = st.d, e = st.e, f = st.f,
  )
end

export Soil, StateBEPS, State2Soil!, state_hydraulic
export VARS_SCALAR, VARS_VECTOR, ALL_VARS_STATE
