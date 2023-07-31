using Biofilm 

# Input parameters
mumax = 20; KM = 3;
p = (
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title =    "Single Solute and Particulate Case",
    tFinal =   1.0,   # Simulation time [days]
    tol =      1e-2,  # Tolerance
    outPeriod =0.1,   # Time between outputs [days]

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames =["Heterotroph"],    # Particulate names
    Xto =   [10.0],        # Tank particulate concentration initial condition(s)
    Pbo =   [0.08],        # Biofilm particulates volume fraction initial condition(s) 
    rho =   [2.0E4],       # Particulate densities
    Kdet =  20000.0,       # Particulates detachment coefficient
    srcX =  [(S,X,Lf,t,z,p) -> 0.0],     # Source of particulates
    # Growthrates for each particulate
    mu = [(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1])],

    # ----------------- #
    # Solute Parameters #
    # ----------------- #
    SNames =["Nutrient"],     # Solute names
    Sin =   [(t) -> 100],   # Solute inflow (can be function of time)
    Sto =   [10.0],         # Tank solute concentration initial condition(s)
    Sbo =   [0.0],          # Biofilm solutes concentration initial condition(s)
    Yxs =   [2.646],        # Biomass yield coefficient on solute
    Dt =    [4.0E-5],       # Aquious solute diffusion through tank fluid
    Db =    [6.9E-5],       # Effective solute diffusion through biofilm
    srcS =  [(S,X,Lf,t,z,p) -> 0.0],     # Source of solutes

    # --------------- #
    # Tank Parameters #
    # --------------- #
    V = 0.1,        # Volume of tank [m³]
    A =   1,        # Surface area of biofilm [m²]
    Q =   1,        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz =  50,       # Number of grid points in biofilm
    Lfo = 1.0E-5,   # Biofilm initial thickness [m]
    LL =  1.00E-7,  # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
biofilm_plot(sol,p) # Plot final results
