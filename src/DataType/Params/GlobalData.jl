export find_SoilType, find_VegType
export SoilTypes, VegTypes
export SOIL_PARAMS

using JSON
const PATH_VEG = joinpath(@__DIR__, "data", "ParamVeg.json")
const PATH_GEN = joinpath(@__DIR__, "data", "ParamGeneral.json")

VegTypes = ["ENF", "DNF", "DBF", "EBF", "Shrub-SH", "CRO", "C4", "default"]
VegCodes = [1, 2, 6, 9, 13, 25, 40, -1]

readGeneralParam() = JSON.parsefile(PATH_GEN)

function veg_name(lc::Integer)
  type_idx = findfirst(==(Int(lc)), VegCodes)
  type_idx === nothing ? "default" : VegTypes[type_idx]
end

veg_param_key(name::AbstractString) = name == "CRO" ? "default" : String(name)

function readVegRaw(lc::Integer=1)
  veg_data = JSON.parsefile(PATH_VEG)
  return veg_data[veg_param_key(veg_name(lc))]
end


const SOIL_PARAMS = [
  (name="sand",
    b=[1.7, 1.9, 2.1, 2.3, 2.5],
    K_sat=[20.88, 18.72, 16.56, 12.6, 3.6],
    θ_sat=0.437, θ_vfc=0.09, θ_res=0.03,
    ψ_sat=[0.07, 0.08, 0.09, 0.10, 0.12],
    κ_dry=8.6),
  (name="loamy_sand",
    b=[2.1, 2.3, 2.5, 2.7, 2.9],
    K_sat=[6.12, 5.4, 5.04, 3.6, 1.08],
    θ_sat=0.437, θ_vfc=0.21, θ_res=0.06,
    ψ_sat=[0.09, 0.10, 0.11, 0.12, 0.14],
    κ_dry=8.3),
  (name="sandy_loam",
    b=[3.1, 3.3, 3.5, 3.7, 3.9],
    K_sat=[2.592, 2.3328, 2.0736, 1.5552, 0.5184],
    θ_sat=0.453, θ_vfc=0.21, θ_res=0.10,
    ψ_sat=[0.15, 0.16, 0.17, 0.18, 0.20],
    κ_dry=8.0),
  (name="loam",
    b=[4.5, 4.7, 4.9, 5.1, 5.3],
    K_sat=[1.332, 1.188, 1.0656, 0.7992, 0.2664],
    θ_sat=0.463, θ_vfc=0.27, θ_res=0.12,
    ψ_sat=[0.11, 0.12, 0.13, 0.14, 0.16],
    κ_dry=7.0),
  (name="silty_loam",
    b=[4.7, 4.9, 5.1, 5.3, 5.5],
    K_sat=[0.684, 0.612, 0.5472, 0.4104, 0.1368],
    θ_sat=0.501, θ_vfc=0.33, θ_res=0.13,
    ψ_sat=[0.21, 0.22, 0.23, 0.24, 0.26],
    κ_dry=6.3),
  (name="sandy_clay_loam",
    b=[4.0, 4.2, 4.4, 4.6, 4.8],
    K_sat=[0.432, 0.3888, 0.3456, 0.2592, 0.0864],
    θ_sat=0.398, θ_vfc=0.26, θ_res=0.15,
    ψ_sat=[0.28, 0.29, 0.30, 0.31, 0.33],
    κ_dry=7.0),
  (name="clay_loam", # 7
    b=[5.2, 5.4, 5.6, 5.8, 6.0],
    K_sat=[0.2304, 0.2088, 0.1836, 0.1368, 0.0468],
    θ_sat=0.464, θ_vfc=0.32, θ_res=0.20,
    ψ_sat=[0.26, 0.27, 0.28, 0.29, 0.31],
    κ_dry=5.8),
  (name="silty_clay_loam",
    b=[6.6, 6.8, 7.0, 7.2, 7.4],
    K_sat=[0.1512, 0.1368, 0.1224, 0.09072, 0.03024],
    θ_sat=0.471, θ_vfc=0.37, θ_res=0.32,
    ψ_sat=[0.33, 0.34, 0.35, 0.36, 0.38],
    κ_dry=4.2),
  (name="sandy_clay",
    b=[6.0, 6.2, 6.4, 6.6, 6.8],
    K_sat=[0.1188, 0.108, 0.09504, 0.07128, 0.02376],
    θ_sat=0.430, θ_vfc=0.34, θ_res=0.24,
    ψ_sat=[0.29, 0.30, 0.31, 0.32, 0.34],
    κ_dry=6.3),
  (name="silty_clay",
    b=[7.9, 8.1, 8.3, 8.5, 8.7],
    K_sat=[0.09, 0.081, 0.072, 0.054, 0.018],
    θ_sat=0.479, θ_vfc=0.39, θ_res=0.25,
    ψ_sat=[0.34, 0.35, 0.36, 0.37, 0.39],
    κ_dry=4.0),
  (name="clay",
    b=[7.6, 7.8, 8.0, 8.2, 8.4],
    K_sat=[0.0612, 0.05508, 0.04896, 0.03672, 0.01224],
    θ_sat=0.475, θ_vfc=0.40, θ_res=0.27,
    ψ_sat=[0.37, 0.38, 0.39, 0.40, 0.42],
    κ_dry=4.4)
]


SoilTypes = [soil.name for soil in SOIL_PARAMS]

function soil_index(code::Integer)
  i = Int(code)
  1 <= i <= length(SOIL_PARAMS) ? i : length(SOIL_PARAMS)
end

soil_name(code::Integer) = SOIL_PARAMS[soil_index(code)].name

find_SoilType(name::AbstractString) = findfirst(SoilTypes .== name)
find_SoilType(code::Integer) = soil_index(code)

find_VegType(code::Integer) = Int(code)

function find_VegType(name::AbstractString)
  I = findfirst(VegTypes .== name)
  isnothing(I) ? -1 : VegCodes[I]
end
