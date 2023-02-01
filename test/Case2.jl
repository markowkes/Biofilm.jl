using Biofilm

# Constants used for growthrates of particulate(s)
mumax = 20;
KM = 3;

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Multiple Independent Substrate Case",
    tFinal=1,       # Simulation time [days]
    tol=1e-2,        # Tolerance
    outPeriod=0.1,   # Time between outputs [days]
    makePlots=false,

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["Bug"],     # Particulate names
    Xto=[10.0],          # Tank particulate concentration initial condition(s)
    Pbo=[0.2],          # Biofilm particulates volume fraction initial condition(s) 
    rho=[1.0e5],        # Particulate densities
    Kdet=1900.0,        # Particulates detachment coefficient
    srcX=[(S,X,t,p) -> 0.0],     # Source of particulates
    # Growthrates for each particulate (constants defined above!)
    mu=[(S,X,Lf,t,z,p) -> (mumax * S[1,:]) ./ (KM .+ S[1,:])],

    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Substrate 1","Substrate 2"],  # Substrate names
    Sin=[(t) -> 25, (t) -> 25], # Substrate inflow (can be function of time)
    Sto=[25.0,25.0],         # Tank substrate concentration initial condition(s)
    Sbo=[0.0,0.0],          # Biofilm substrates concentration initial condition(s)
    Yxs=[0.5 0.0],          # Biomass yield coefficient on substrate
    Daq=[4.0e-5, 6.0e-5],   # Substrate diffusion through boundary layer
    De =[1.0e-5, 1.5e-5],   # Substrate diffusion through biofilm     
    srcS=[(S,X,t,p) -> 0.0,
          (S,X,t,p) -> 0.0],     # Source of substrates
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.1,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=1,          # Flowrate through tank [m³/s]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=50,          # Number of grid points in biofilm
    Lfo=5.0e-6,     # Biofilm initial thickness [m]
    LL=1.0e-4,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver