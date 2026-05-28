# JAX 风格 setup：st, ps = setup(model)
# st = 状态变量(会变), ps = 模型参数(不变)
export setup


## Internal Helpers ============================================================
# 从 ParamBEPS 初始化 StateBEPS
function _calc_Δz₊ₕ_cm(dz_cm::AbstractVector, N::Int)
  Δz₊ₕ = zeros(10)
  for i in 1:N-1
    Δz₊ₕ[i] = 0.5 * (dz_cm[i] + dz_cm[i+1])
  end
  N > 0 && (Δz₊ₕ[N] = 0.5 * dz_cm[N])
  Δz₊ₕ
end

function _init_state(ps::ParamBEPS, Tsoil, Ta, θ0, z_snow)
  N = ps.N
  st = StateBEPS(; n_layer=Cint(N), N)
  st.dz[1:N] .= ps.dz
  st.Δz_cm[1:N] .= ps.hydraulic.dz_cm[1:N]
  st.Δz₊ₕ_cm   .= _calc_Δz₊ₕ_cm(ps.hydraulic.dz_cm, N)
  UpdateRootFraction!(st, ps)
  Init_Soil_T_θ!(st, Float64(Tsoil), Float64(Ta), Float64(θ0), Float64(z_snow))
  st.Tsnow_c .= Float64(Ta)
  st
end

# 初始化 Soil (Julia/C 版本)
function _init_soil(version::String, VegType::Union{AbstractString,Integer}, SoilType::Union{AbstractString,Integer},
    r_drainage, r_root_decay, Tsoil, Ta, θ0, z_snow)
  soil = version == "julia" ? Soil() : Soil_c()
  soil.r_drainage = r_drainage

  VegCode = find_VegType(VegType)
  SoilCode = find_SoilType(SoilType)
  
  Init_Soil_Parameters(soil, VegCode, SoilCode, r_root_decay)
  Init_Soil_T_θ!(soil, Float64(Tsoil), Float64(Ta), Float64(θ0), Float64(z_snow))
  soil
end


## Setup Functions (JAX Style) =================================================
"从 Soil 拆分为 (StateBEPS, ParamBEPS)"
function setup(soil::Soil; FT=Float64)
  st = StateBEPS(soil)
  ps = ParamBEPS{FT}(); Soil2Params!(ps, soil)
  st, ps
end

"直接构造 (StateBEPS, ParamBEPS)，无需先创建 Soil"
function setup(VegType::Union{AbstractString,Integer}, SoilType::Union{AbstractString,Integer};
    Ta::Real=20.0, Tsoil::Real=Ta, θ0::Real=0.3, z_snow::Real=0.0,
    r_drainage::Real=0.5, r_root_decay::Real=0.95, N::Int=5, FT::Type=Float64)
  ps = ParamBEPS(VegType, SoilType; N, FT, r_drainage)
  ps.veg.r_root_decay = FT(r_root_decay)
  _init_state(ps, Tsoil, Ta, θ0, z_snow), ps
end

"纯关键字参数版本"
setup(; VegType, SoilType, kw...) = setup(VegType, SoilType; kw...)

"从已有 ParamBEPS 构造状态"
function setup(ps::ParamBEPS; Ta::Real=20.0, Tsoil::Real=Ta, θ0::Real=0.3, z_snow::Real=0.0)
  _init_state(ps, Tsoil, Ta, θ0, z_snow), ps
end


## Legacy Setup (setup_model) ==================================================
"统一 setup，通过 version 区分 Julia/C 版本，返回 (Soil, State, Params)"
function setup_model(VegType::Union{AbstractString,Integer}, SoilType::Union{AbstractString,Integer};
    version::String="julia", Ta::Real=20.0, Tsoil::Real=Ta, θ0::Real=0.3, z_snow::Real=0.0,
    r_drainage::Real=0.5, r_root_decay::Real=0.95, FT::Type=Float64)
  soil = _init_soil(version, VegType, SoilType, r_drainage, r_root_decay, Tsoil, Ta, θ0, z_snow)
  state = version == "julia" ? StateBEPS(soil) : zeros(41)
  InitState!(soil, state, Float64(Ta))
  ps = ParamBEPS(VegType, SoilType; FT)
  version == "julia" && Soil2Params!(ps, soil)
  soil, state, ps
end

setup_jl(args...; kw...) = setup_model(args...; version="julia", kw...)
setup_c(args...; kw...) = setup_model(args...; version="c", kw...)

export setup_model, setup_c, setup_jl
