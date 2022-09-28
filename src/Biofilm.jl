module Biofilm

export param, BiofilmSolver, makePlots, analyzeBiofilm, movieBiofilm

include("outputs.jl")
include("structs.jl")
include("rhs.jl")
include("tools.jl")
include("solver.jl")
include("checkParameters.jl")
include("computes.jl")

end