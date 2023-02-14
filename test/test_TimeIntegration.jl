# test_TimeIntegration.jl
# -----------------------------------------------
# Checks the temporal changes in substrate in 
# tank without a particulate matches analytical solution
# -----------------------------------------------

using Biofilm
using Statistics
using Printf
#using Plots

function test_TimeIntegration()

    # Constants used for growthrates of particulate(s)
    mumax = 2000;
    KM = 2500;

    # Define a structure to hold all the parameters
    p = param(
        # --------------------- #
        # Simulation Parameters #
        # --------------------- #
        Title="Time Integration Test",
        tFinal=1,       # Simulation time [days]
        tol=1e-8,       # Tolerance
        outPeriod=1.0,  # Time between outputs [days]
        makePlots=false,

        # ---------------------- #
        # Particulate Parameters #
        # ---------------------- #
        XNames=["Bug"],     # Particulate names
        Xto=[10.0],          # Tank particulate concentraion initial condition(s)
        Pbo=[0.08],         # Biofilm particulates volume fraction initial condition(s) 
        rho=[2.0E4],        # Particulate densities
        Kdet=20000.0,       # Particulates detachment coefficient
        srcX=[(S,X,Lf,t,z,p) -> 0.0],     # Source of particulates
        # Growthrates for each particulate (constants defined above!)
        mu=[(S,X,Lf,t,z,p) -> 0.0.*S[1]],

        # -------------------- #
        # Substrate Parameters #
        # -------------------- #
        SNames=["Oxygen"],   # Substrate names
        Sin=[(t) -> 100],    # Substrate inflow (can be function of time)
        Sto=[25.0],          # Tank substrate concentraion initial condition(s)
        Sbo=[0.0],           # Biofilm substrates concentration initial condition(s)
        Yxs=[0.5],           # Biomass yield coefficient on substrate
        Daq=[1.0E-15],        # Substrate diffusion through boundary layer
        De =[1.0E-15],        # Substrate diffusion through biofilm     
        srcS=[(S,X,Lf,t,z,p) -> 0.0],     # Source of substrates
        
        # --------------- #
        # Tank Parameters #
        # --------------- #
        V=0.1,        # Volume of tank [m³]
        A=1,          # Surface area of biofilm [m²]
        Q=1,        # Flowrate through tank [m³/d]

        # ------------------ #
        # Biofilm Parameters #
        # ------------------ #
        Nz=20,          # Number of grid points in biofilm
        Lfo=5.0E-5,     # Biofilm initial thickness [m]
        LL=0.0,         # Boundary layer thickness [m]
    )

    t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
    
    # Analytic Result
    Sin = p.Sin[1](t[end])
    St_ana = Sin .+ (p.Sto[1] - Sin).*exp.(-p.Q/p.V.*t)

    # Plots of solutions and error
    # fig = plot(t,St',label="Simulated")
    # fig = plot!(t,St_ana,label="Analytic",line=:dash)
    # display(fig)

    # fig = plot(t,abs.(St'.-St_ana),label="Error")
    # display(fig)

    error = maximum(abs.(St_ana .- St'))

    println("Maximum relative error = ",error/Sin)

    return error
end

error = test_TimeIntegration()