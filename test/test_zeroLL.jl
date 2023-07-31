# test_zeroLL.jl
# -----------------------------------------------
# Checks the ability of the code to run when LL=0.
# Expected output is that S_top = St (tank)
# -----------------------------------------------

using Biofilm

function test_zeroLL()

    # Constants used for growthrates of particulate(s)
    mumax = 20;
    KM = 3;

    # Define a structure to hold all the parameters
    p = (
        # --------------------- #
        # Simulation Parameters #
        # --------------------- #
        Title="Zero LL Test",
        tFinal=1,       # Simulation time [days]
        tol=1e-2,       # Tolerance
        outPeriod=0.1,  # Time between outputs [days]
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
        mu=[(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1])],

        # -------------------- #
        # Solute Parameters    #
        # -------------------- #
        SNames=["Oxygen"],  # Solute names
        Sin=[(t) -> 100],   # Solute inflow (can be function of time)
        Sto=[10.0],          # Tank solute concentraion initial condition(s)
        Sbo=[0.0],          # Biofilm solutes concentration initial condition(s)
        Yxs=[2.646],        # Biomass yield coefficient on solute
        Dt=[4.0E-5],       # Aquious solute diffusion through tank fluid
        Db=[6.9E-5],        # Effective solute diffusion through biofilm
        srcS=[(S,X,Lf,t,z,p) -> 0.0],     # Source of solutes
        
        # --------------- #
        # Tank Parameters #
        # --------------- #
        V=0.1,        # Volume of tank [m³]
        A=1,          # Surface area of biofilm [m²]
        Q=1,          # Flowrate through tank [m³/d]

        # ------------------ #
        # Biofilm Parameters #
        # ------------------ #
        Nz=50,          # Number of grid points in biofilm
        Lfo=1.0E-5,     # Biofilm initial thickness [m]
        LL=0.0,   # Boundary layer thickness [m]
    )

    t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver

    # Compute S_top
    g = computeGrid(Lf[end],p)
    S_top = Biofilm.computeS_top(St[:,end],Sb,p,g)

    # Final values
    computed = S_top[1]
    analytic = St[end]

    return computed, analytic
end

computed, analytic = test_zeroLL()