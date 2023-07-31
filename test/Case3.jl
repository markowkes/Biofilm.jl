using Biofilm 

# Input constants 
b = 0.1 # Source term constant
mumax = 2; KM = 1; # Growthrate constants

# Tuple to hold all input parameters
p = (
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title =   "Multiple Particulate Case",
    tFinal =  100,   # Simulation time [days]
    tol =     1e-4,  # Tolerance
    outPeriod = 5.0, # Time between outputs [days]
    makePlots = false,

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames = ["Living Bug","Dead Bug"],    # Particulate names
    Xto =  [10.0, 0.0],        # Tank particulate concentration initial condition(s)
    Pbo =  [0.08, 0.0],        # Biofilm particulates volume fraction initial condition(s) 
    rho =  [2.0e5, 2.0e5],     # Particulate densities
    Kdet = 1980.0,             # Particulates detachment coefficient
    srcX = [(S,X,Lf,t,z,p) -> -b*X[1], # Source of particulates
            (S,X,Lf,t,z,p) ->  b*X[1]],
    # Growthrates for each particulate
    mu =[(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1]) 
         (S,X,Lf,t,z,p) -> 0.0 ],

    # ----------------- #
    # Solute Parameters #
    # ----------------- #
    SNames = ["Solute"], # Solute names
    Sin =  [(t) -> 25],     # Solute inflow (can be function of time)
    Sto =  [25.0],          # Tank solute concentration initial condition(s)
    Sbo =  [0.0],           # Biofilm solutes concentration initial condition(s)
    Yxs =  [0.378, 0],      # Biomass yield coefficient on solute
    Dt =   [1.38e-4],       # Aquious solute diffusion through tank fluid
    Db =   [6.9E-5],        # Effective solute diffusion through biofilm
    srcS = [(S,X,Lf,t,z,p) -> 0.0], # Source of solutes

    # --------------- #
    # Tank Parameters #
    # --------------- #
    V =0.1,        # Volume of tank [m³]
    A =  1,        # Surface area of biofilm [m²]
    Q =  1,        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz = 50,       # Number of grid points in biofilm
    Lfo =5.0e-6,   # Biofilm initial thickness [m]
    LL = 1.0e-5,  # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver