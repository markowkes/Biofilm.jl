module Biofilm

export BiofilmSolver, computeGrid
export biofilm_analyze, biofilm_movie, MeanBiofilmVarsWithTime, biofilm_sol2csv

include("outputs.jl")
include("structs.jl")
include("rhs.jl")
include("tools.jl")
include("solver.jl")
include("parameters.jl")
include("computes.jl")
include("postprocess.jl")
include("precompile.jl")

function julia_main()::Cint

    # Load and process input file
    p = process_command_line_args(ARGS)

    # Run solver
    t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) 

    return 0 # if things finished successfully
end

end