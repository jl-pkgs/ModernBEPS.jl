@bounds @with_kw_noshow mutable struct ParamBEPS{FT<:AbstractFloat}
  N::Int = 5
  dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]  # 土壤层厚度 [m], BEPS V2023
  r_drainage::FT = Cdouble(0.50) | (0.2, 0.7)  # ? 地表排水速率（地表汇流），可考虑采用曼宁公式

  ψ_min::FT = Cdouble(33.0)  # [m], about 0.10~0.33 MPa开始胁迫点
  alpha::FT = Cdouble(0.4)   # [-], 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  hydraulic::ParamSoilHydraulicLayers{FT} = ParamSoilHydraulicLayers{FT,N}()
  thermal::ParamSoilThermalLayers{FT} = ParamSoilThermalLayers{FT,N}()

  veg::ParamVeg{FT} = ParamVeg{FT}()
end

# `kw...`: other params like, `r_drainage`
function ParamBEPS(VegType::Union{AbstractString,Integer}, SoilType::Union{AbstractString,Integer};
  N::Int=5, FT=Float64, kw...)

  veg = InitParam_Veg(VegType; FT)
  hydraulic, thermal = InitParam_Soil(SoilType, N, FT)

  ψ_min = veg.is_bforest ? FT(10.0) : FT(33.0) # 开始胁迫点
  alpha = veg.is_bforest ? FT(1.5) : FT(0.4)   # 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  ParamBEPS{FT}(;
    N, kw..., ψ_min, alpha,
    hydraulic, thermal, veg
  )
end


# 这里应该加一个show function，打印模型参数信息
function Base.show(io::IO, model::M) where {M<:ParamBEPS}
  printstyled(io, "$M, N = $(model.N)\n", color=:blue, bold=true)

  fields_all = fieldnames(M)
  fields = setdiff(fields_all, [:N, :hydraulic, :thermal, :veg])

  n = length(fields)
  for i = 1:n
    field = fields[i]
    value = getfield(model, field)
    type = typeof(value)
    isa(value, Function) && (type = Function)
    println(io, "  $field\t: {$type} $value")
    # (i != n) && print(io, "\n")
  end

  ss = 60
  println(io, "-"^ss)
  printstyled(io, "Hydraulic: ", color=:blue, bold=true)
  print(io, model.hydraulic)

  println(io, "-"^ss)
  printstyled(io, "Thermal: ", color=:blue, bold=true)
  print(io, model.thermal)

  println(io, "-"^ss)
  printstyled(io, "Veg: ", color=:blue, bold=true)
  print(io, model.veg)
  print("-"^ss)
  return nothing
end


# DBF or EBF, low constaint threshold
function Params2Soil!(soil::Soil, params::ParamBEPS{FT}; BF=false) where {FT}
  soil.ψ_min = BF ? 10.0 : 33.0 # [m], about 0.10~0.33 MPa开始胁迫点
  soil.alpha = BF ? 1.5 : 0.4   # He 2017 JGR-B, Eq. 4

  (; hydraulic, thermal, N) = params
  soil.n_layer = Cint(N)
  soil.dz[1:5] .= [0.05, 0.10, 0.20, 0.40, 1.25] # BEPS V2023, 土壤层厚度[m]

  soil.r_drainage = Cdouble(params.r_drainage)
  soil.r_root_decay = Cdouble(params.veg.r_root_decay)
  UpdateRootFraction!(soil) # 更新根系分布

  soil.ψ_min = Cdouble(params.ψ_min)
  soil.alpha = Cdouble(params.alpha)

  # soil.θ_vfc[1:N] .= Cdouble.(hydraulic.θ_vfc)
  soil.θ_vwp[1:N] .= Cdouble.(hydraulic.θ_vwp)
  soil.θ_sat[1:N] .= Cdouble.(hydraulic.θ_sat)
  soil.K_sat[1:N] .= Cdouble.(hydraulic.K_sat)
  soil.ψ_sat[1:N] .= Cdouble.(hydraulic.ψ_sat)
  soil.b[1:N] .= Cdouble.(hydraulic.b)

  soil.κ_dry[1:N] .= Cdouble.(thermal.κ_dry)
  soil.ρ_soil[1:N] .= Cdouble.(thermal.ρ_soil)
  soil.V_SOM[1:N] .= Cdouble.(thermal.V_SOM)
end
Params2Soil!(soil::AbstractSoil, params::Nothing) = nothing


function Soil2Params!(params::ParamBEPS{FT}, soil::Soil) where {FT}
  N = Int(soil.n_layer)
  params.N = N

  params.r_drainage = FT(soil.r_drainage)
  params.veg.r_root_decay = FT(soil.r_root_decay)
  params.ψ_min = FT(soil.ψ_min)
  params.alpha = FT(soil.alpha)

  if length(params.dz) != N
    resize!(params.dz, N)
  end
  params.dz .= FT.(soil.dz[1:N])

  (; hydraulic, thermal) = params

  # :θ_vfc
  for field in (:θ_vwp, :θ_sat, :K_sat, :ψ_sat, :b)
    dest = getfield(hydraulic, field)
    src = getfield(soil, field)
    if length(dest) != N
      resize!(dest, N)
    end
    dest .= FT.(src[1:N])
  end

  for field in (:κ_dry, :ρ_soil, :V_SOM)
    dest = getfield(thermal, field)
    src = getfield(soil, field)
    if length(dest) != N
      resize!(dest, N)
    end
    dest .= FT.(src[1:N])
  end
  return params
end




using Ipaper: approx

function interp_depths(SM::AbstractMatrix{FT}, z_obs) where {FT<:AbstractFloat}
  dz::Vector{FT} = FT[0.05, 0.10, 0.20, 0.40, 1.25]
  z_sim = cumsum(dz) .- dz ./ 2

  ntime = size(SM, 1)
  R = zeros(FT, ntime, length(z_obs))
  for i in 1:ntime
    y_sim = @view SM[i, :]
    R[i, :] .= approx(z_sim, y_sim, z_obs)
  end
  R
end

export interp_depths
