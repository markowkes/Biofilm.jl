module Biofilm

export param, BiofilmSolver, computeGrid
export makePlots, analyzeBiofilm, movieBiofilm, MeanBiofilmVarsWithTime, sol2csv
export createDict, addParam!, packageCheckParam, printDict

include("outputs.jl")
include("structs.jl")
include("rhs.jl")
include("tools.jl")
include("solver.jl")
include("parameters.jl")
include("computes.jl")
include("postprocess.jl")

end