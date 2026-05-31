module BEPS
# using BEPS

using DocStringExtensions
using UnPack
import Parameters: @with_kw, @with_kw_noshow
using Printf
using Reexport
using Dates
import DataFrames: DataFrame

using StaticArrays
using Statistics
using ModelParams
using ModelParams: of_RMSE

@reexport using Serialization: deserialize, serialize
@reexport using DelimitedFiles: readdlm
export beps_main
export beps_modern
export Soil_c

export path_proj
path_proj(f...) = normpath(joinpath(@__DIR__, "..", f...))

using LazyArtifacts
using Libdl: dlext
using Base.BinaryPlatforms: HostPlatform, os, arch

const libbeps = let
  _ext  = dlext
  _os   = os(HostPlatform())
  _ar   = arch(HostPlatform()) == "aarch64" ? "arm64" : "x86_64"

  fname = "libbeps-$_os-$_ar.$_ext"
  local_lib = path_proj("deps", fname)
  isfile(local_lib) ? local_lib : joinpath(artifact"libbeps", fname)
end

# import Statistics: mean, std
# include("DataFrames.jl")
# include("Ipaper.jl")
# include("c2julia.jl")
include("SPAC/SPAC.jl")
include("DataType/DataType.jl")
include("BEPS_modules.jl")

include("clang/BEPS_c.jl")
@reexport import BEPS.clang;
import BEPS.clang: inter_prg_c, photosynthesis_c, Soil_c,
  snowpack_stage1, snowpack_stage2, snowpack_stage3
using BEPS.clang

# include("beps_inter_prg.jl")
include("beps_main.jl")
include("beps_modern.jl")
include("ultilize.jl")
include("SoilWater.jl")
include("Optim.jl")

end # module BEPS
