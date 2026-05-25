# 快为第一准则，先把模型跑起来
using BEPS

# 初始化模型参数
VegType::Int = 25
SoilType::Int = 8

model = ParamBEPS(VegType, SoilType)

## 初始化状态变量
Ta = 20.0
Tsoil0 = Ta
θ0 = model.hydraulic.θ_sat[1] .* 0.6
z_snow0 = 0.0

state, model = setup(model; Ta, Tsoil=Float64(Tsoil0), θ0=Float64(θ0), z_snow=Float64(z_snow0))
state

## 开干
df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates; ps=model, state)
