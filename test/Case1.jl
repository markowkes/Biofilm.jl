using Biofilm

# Constants used for growthrates of particulate(s)
mumax = 20;
KM = 3;

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Single Substrate and Particulate Case",
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
    mu=[(S,X,Lf,t,z,p) -> (mumax * S[1,:]) ./ (KM .+ S[1,:])],

    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Oxygen"],  # Substrate names
    Sin=[(t) -> 100],   # Substrate inflow (can be function of time)
    Sto=[10.0],          # Tank substrate concentraion initial condition(s)
    Sbo=[0.0],          # Biofilm substrates concentration initial condition(s)
    Yxs=[2.646],        # Biomass yield coefficient on substrate
    Daq=[4.0E-5],       # Substrate diffusion through boundary layer
    De=[6.9E-5],        # Substrate diffusion through biofilm     
    srcS=[(S,X,Lf,t,z,p) -> 0.0],     # Source of substrates
    
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
    LL=1.00E-7,   # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver