using BEPS, Test

include("test-beps_main.jl")

include("test-macro.jl")
include("test-StateSeries.jl")
include("dev/test-aerodynamic_conductance.jl")
include("test-photosynthesis_standalone.jl")

include("test-beps_modern.jl")
include("test-soilwater.jl")

include("test-ModelParams.jl")
include("test-utilize.jl")
include("test-UpdateSoilMoisture_Q0.jl")
include("modules/modules.jl")
