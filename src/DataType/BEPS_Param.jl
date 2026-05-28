const _BEPS_DZ = (0.05, 0.10, 0.20, 0.40, 1.25)

function _default_dz(::Type{FT}, N::Int) where {FT<:AbstractFloat}
  base = FT.(_BEPS_DZ)
  N <= length(base) && return collect(base[1:N])
  return [collect(base); fill(last(base), N - length(base))]
end


@bounds @with_kw_noshow mutable struct ParamBEPS{FT<:AbstractFloat,N,
  H<:HydraulicProfile{FT,N},T<:ThermalProfile{FT,N}} <: AbstractSoilModel{FT,N}
  dz::Vector{FT} = _default_dz(FT, N)
  r_drainage::FT = FT(0.50) | (0.2, 0.7)

  ψ_min::FT = FT(33.0)
  alpha::FT = FT(0.4)

  hydraulic::H
  thermal::T
  veg::ParamVeg{FT} = ParamVeg{FT}()
end


function ParamBEPS{FT,N}(hydraulic::H, thermal::T; dz=_default_dz(FT, N), kwargs...) where {
  FT<:AbstractFloat,N,H<:HydraulicProfile{FT,N},T<:ThermalProfile{FT,N}}

  dz = FT.(collect(dz))
  hydraulic.dz_cm = FT.(100 .* dz)
  ParamBEPS{FT,N,H,T}(; dz, hydraulic, thermal, kwargs...)
end

function ParamBEPS{FT,N}(; dz=_default_dz(FT, N), kwargs...) where {FT<:AbstractFloat,N}
  dz = FT.(collect(dz))
  hydraulic = HydraulicProfile{FT,N}(CampbellLayers{FT,N}(), KvLayers{FT,N}(), FT.(100 .* dz))
  thermal = ThermalProfile{FT,N}(ThermalBaseLayers{FT,N}())
  ParamBEPS{FT,N}(hydraulic, thermal; dz, kwargs...)
end

ParamBEPS{FT}(; N::Int=5, kwargs...) where {FT<:AbstractFloat} = ParamBEPS{FT,N}(; kwargs...)

# `kw...`: other params like, `r_drainage`
function ParamBEPS(VegType::Union{AbstractString,Integer}, SoilType::Union{AbstractString,Integer};
  N::Int=5, FT::Type{<:AbstractFloat}=Float64, kw...)

  veg = InitParam_Veg(VegType; FT)
  hydraulic, thermal = InitParam_Soil(SoilType, N, FT)

  ψ_min = veg.is_bforest ? FT(10.0) : FT(33.0) # 开始胁迫点
  alpha = veg.is_bforest ? FT(1.5) : FT(0.4)   # 土壤水限制因子参数，He 2017 JGR-B, Eq. 4

  ParamBEPS{FT,N}(hydraulic, thermal;
    kw..., ψ_min, alpha, veg
  )
end

nlayer(::ParamBEPS{FT,N}) where {FT,N} = N

function Base.getproperty(x::ParamBEPS, name::Symbol)
  name === :N && return nlayer(x)
  name === :r_root_decay && return getfield(x, :veg).r_root_decay
  return getfield(x, name)
end


_fit_layers(::Type{FT}, values, N::Int) where {FT<:AbstractFloat} = begin
  xs = FT.(collect(values))
  length(xs) >= N ? xs[1:N] : [xs; fill(last(xs), N - length(xs))]
end


function Base.getproperty(x::HydraulicProfile, name::Symbol)
  name in (:profile, :layers, :kv, :dz_cm) && return getfield(x, name)
  profile = getfield(x, :profile)
  name === :K_sat && return getproperty(profile, :Ksat)
  return getproperty(profile, name)
end

function Base.getproperty(x::ThermalProfile, name::Symbol)
  name in (:profile, :layers) && return getfield(x, name)
  return getproperty(getfield(x, :profile), name)
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
  soil.dz[1:N] .= params.dz[1:N]

  soil.r_drainage = Cdouble(params.r_drainage)
  soil.r_root_decay = Cdouble(params.veg.r_root_decay)
  UpdateRootFraction!(soil) # 更新根系分布

  soil.ψ_min = Cdouble(params.ψ_min)
  soil.alpha = Cdouble(params.alpha)

  # soil.θ_vfc[1:N] .= Cdouble.(hydraulic.θ_vfc)
  soil.θ_res[1:N] .= Cdouble.(hydraulic.θ_res)
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
  N == params.N || throw(ArgumentError("Cannot write $N soil layers into ParamBEPS{$FT,$(params.N)}."))

  params.r_drainage = FT(soil.r_drainage)
  params.veg.r_root_decay = FT(soil.r_root_decay)
  params.ψ_min = FT(soil.ψ_min)
  params.alpha = FT(soil.alpha)

  if length(params.dz) != N
    resize!(params.dz, N)
  end
  params.dz .= FT.(soil.dz[1:N])

  (; hydraulic, thermal) = params

  hydraulic.profile.θ_res .= FT.(soil.θ_res[1:N])
  hydraulic.profile.θ_sat .= FT.(soil.θ_sat[1:N])
  hydraulic.profile.Ksat .= FT.(soil.K_sat[1:N])
  hydraulic.kv.kv .= FT.(soil.K_sat[1:N])
  hydraulic.profile.ψ_sat .= FT.(soil.ψ_sat[1:N])
  hydraulic.profile.b .= FT.(soil.b[1:N])
  hydraulic.dz_cm .= FT.(100 .* params.dz)

  thermal.profile.κ_dry .= FT.(soil.κ_dry[1:N])
  thermal.profile.ρ_soil .= FT.(soil.ρ_soil[1:N])
  thermal.profile.V_SOM .= FT.(soil.V_SOM[1:N])

  for i in 1:N
    hydraulic.layers[i] = hydraulic.profile[i]
    thermal.layers[i] = thermal.profile[i]
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
