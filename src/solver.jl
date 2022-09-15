# ODE Progress
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using DifferentialEquations
using UnPack

function BiofilmSolver(p::param)

    println("Starting Solver ...")

    # Check parameters passed to solver 
    checkParameters(p)

    # Unpack parameters
    @unpack Nx,Ns,Nz,Xo,So,Pbo,Sbo,Lfo,tol,tFinal,outPeriod,discontinuityPeriod = p

    # Compute ranges of dependent variables in sol array
    # sol=[X,S,Pb,S,Lf]
    nVar=0; 
    N=Nx;    rX =nVar+1:nVar+N; nVar+=N # X =u[rX]
    N=Ns;    rS =nVar+1:nVar+N; nVar+=N # S =u[rS]
    N=Nx*Nz; rPb=nVar+1:nVar+N; nVar+=N # Pb=u[rPb]
    N=Ns*Nz; rSb=nVar+1:nVar+N; nVar+=N # Sb=u[rSb]
    N=1;     rLf=nVar+1:nVar+N          # Lf=u[rLf]

    # Store ranges in struct to be used in RHS calc
    r = ranges(rX,rS,rPb,rSb,rLf)

    # Prepare biofilm initial conditions
    Pbo_grid = vec(Pbo'.*ones(Nz,Nx))
    Sbo_grid = vec(Sbo'.*ones(Nz,Ns))
    
    # Package initial conditions as one continuous vector
    sol0=vcat(Xo,So,Pbo_grid,Sbo_grid,Lfo)

    # Prepare ODE Solver 
    prob = ODEProblem(biofilmRHS!,sol0,(0.0,p.tFinal),[p,r])

    # Compute solver step size 
    # greatest common divisor of output and discontinutity periods
    solverStep=gcd(outPeriod,discontinuityPeriod)
    
    # Output times 
    outTimes = range(start=0.0,step=solverStep,stop=tFinal)
    affect!(integrator) = outputs(integrator)
    cb = PresetTimeCallback(outTimes,affect!)
    
    # Run solver
    sol=solve(prob,
        reltol=tol,
        abstol=tol,
        callback=cb,
        progress = true,
        progress_steps = 100,
        alg_hints = [:stiff],
        )

    # Convert solution to dependent variables
    t,X,S,Pb,Sb,Lf=unpack_solution(sol,p,r)

    return sol,t,X,S,Pb,Sb,Lf   
    #return sol
end