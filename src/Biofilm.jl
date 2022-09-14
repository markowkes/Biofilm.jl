module Biofilm

export param, BiofilmSolver,outputs

include("outputs.jl")
include("structs.jl")
include("rhs.jl")
include("tools.jl")
include("solver.jl")
include("checkParameters.jl")

end