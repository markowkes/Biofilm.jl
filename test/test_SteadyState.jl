# test_SteadyState.jl
# -----------------------------------------------
# Test with large diffusivity with an analytic solution
# -----------------------------------------------

using Biofilm
using Statistics
using Printf

function test_SteadyState()
    # Constants used for growthrates of particulate(s)
    mumax = 2000;
    KM = 2500;

    # Define a structure to hold all the parameters
    p = (
        # --------------------- #
        # Simulation Parameters #
        # --------------------- #
        Title="Diffusion Test",
        tFinal=10,       # Simulation time [days]
        tol=1e-4,       # Tolerance
        outPeriod=1.0,  # Time between outputs [days]
        makePlots=false,

        # ---------------------- #
        # Particulate Parameters #
        # ---------------------- #
        XNames=["Bug"],     # Particulate names
        Xto=[10.0],          # Tank particulate concentraion initial condition(s)
        Pbo=[0.08],         # Biofilm particulates volume fraction initial condition(s) 
        rho=[2.0E4],        # Particulate densities
        Kdet=1900.0,       # Particulates detachment coefficient
        srcX=[(S,X,Lf,t,z,p) -> 0.0],     # Source of particulates
        # Growthrates for each particulate (constants defined above!)
        mu=[(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM)],

        # -------------------- #
        # Solute Parameters    #
        # -------------------- #
        SNames=["Oxygen"],   # Solute names
        Sin=[(t) -> 100],    # Solute inflow (can be function of time)
        Sto=[25.0],          # Tank solute concentraion initial condition(s)
        Sbo=[0.0],           # Biofilm solutes concentration initial condition(s)
        Yxs=[0.5],           # Biomass yield coefficient on solute
        Dt=[100.0],         # Aquious solute diffusion through tank fluid
        Db=[100.0],         # Effective solute diffusion through biofilm
        srcS=[(S,X,Lf,t,z,p) -> 0.0],     # Source of solutes
        
        # --------------- #
        # Tank Parameters #
        # --------------- #
        V=0.1,        # Volume of tank [m³]
        A=1,          # Surface area of biofilm [m²]
        Q=1,        # Flowrate through tank [m³/d]

        # ------------------ #
        # Biofilm Parameters #
        # ------------------ #
        Nz=50,          # Number of grid points in biofilm
        Lfo=5.0E-6,     # Biofilm initial thickness [m]
        LL=1e-4,         # Boundary layer thickness [m]
    )
    
    t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver

    # Analytic Result
    St_ana=[1.0]
    Xb_ana = mean(Pb)*p.rho[1]
    Res = 1.0
    g = Biofilm.computeGrid(Lf[end],p)
    while abs(Res) > 1e-10
        # Note that only value of St matters due to definition of mu
        mu_ana = p.mu[1](St_ana,Xt[:,end],Lf[end],t[end],zm,p)[1]
        Lf_ana = mu_ana / p.Kdet
        Xt_ana = p.Yxs[1]*p.Q/(mu_ana*p.V)*(p.Sin[1](t[end]) - St_ana[1]) - Xb_ana*Lf_ana*p.A/p.V
        Res = mu_ana*Xt_ana - p.Q*Xt_ana/p.V + mu_ana*Lf_ana*p.A*Xb_ana/p.V 
        St_ana[1] -= 0.001*Res
    end

    # Final values
    computed = St[1,end]
    analytic = St_ana[1]

    # % difference 
    println("% difference = ",abs(computed - analytic)/analytic*100)

    return computed, analytic
end

computed, analytic = test_SteadyState()