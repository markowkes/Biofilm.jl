module Biofilm

export param, BiofilmSolver, computeGrid, makePlots, analyzeBiofilm, movieBiofilm, MeanBiofilmVarsWithTime, sol2csv

include("outputs.jl")
include("structs.jl")
include("rhs.jl")
include("tools.jl")
include("solver.jl")
include("checkParameters.jl")
include("computes.jl")
include("postprocess.jl")

end