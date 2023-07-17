# ODE Progress
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#using DifferentialEquations
using OrdinaryDiffEq
using Sundials
using UnPack
using DiffEqCallbacks

"""
    BiofilmSolver(p::param)

Take parameters defining a biofilm case and computes solution
"""
function BiofilmSolver(p)

    println("Starting Solver ...")

    # Check inputs 
    p = checkInputs(p)

    # Unpack parameters
    @unpack Nx,Ns,Nz,Xto,Sto,Pbo,Sbo,Lfo,tol,tFinal,outPeriod,discontinuityPeriod = p

    # Compute ranges of dependent variables in sol array (add to p)
    p = computeRanges(p)

    # Prepare biofilm initial conditions
    Pbo_grid = vec(Pbo.*ones(Nx,Nz))
    Sbo_grid = vec(Sbo.*ones(Ns,Nz))
    
    # Package initial conditions as one continuous vector
    sol0=vcat(Xto,Sto,Pbo_grid,Sbo_grid,Lfo)

    # Prepare ODE Solver 
    prob = ODEProblem(biofilmRHS!,sol0,(0.0,p.tFinal),p)

    # Compute solver step size 
    # greatest common divisor of output and discontinutity periods
    solverStep=gcd(outPeriod,discontinuityPeriod)
    
    # Output times 
    outTimes = range(start=0.0,step=solverStep,stop=tFinal)
    affect!(integrator) = outputs(integrator)
    cb = PresetTimeCallback(outTimes,affect!)
    
    # Run solver
    #GC.enable(false)
    sol=solve(prob,
        Sundials.CVODE_BDF(),
        reltol=tol,
        abstol=tol,
        callback=cb,
        progress = true,
        progress_steps = 100,
        alg_hints = [:stiff],
        
        )
    #GC.enable(true)

    # Convert solution to dependent variables
    t,Xt,St,Pb,Sb,Lf=unpack_solutionForPlot(sol,p)

    # Final biofilm grid
    z=range(0.0,Lf[end],Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])

    return t,zm,Xt,St,Pb,Sb,Lf,sol
end