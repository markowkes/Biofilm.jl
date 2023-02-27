module Biofilm

export param, BiofilmSolver, computeGrid
export biofilm_analyze, biofilm_movie, MeanBiofilmVarsWithTime, biofilm_sol2csv
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