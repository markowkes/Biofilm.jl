using Biofilm 

# Constants used for growthrates of particulate(s)
mumax = 2;
KM = 1;

# Source term constants
b=0.1

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Multiple Particulate Case",
    tFinal=100,     # Simulation time [days]
    tol=1e-4,       # Tolerance
    outPeriod=5,    # Time between outputs [days]
    optionalPlot="source", # 6th plot: "growthrate" (default) or "source"
    
    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["Living Bug","Dead Bug"], # Particulate names
    Xto=[10.0,0.0],  # Tank particulate concentration initial condition(s)
    Pbo=[0.08,0.0], # Biofilm particulates volume fraction initial condition(s) 
    rho=[2.5e5,2.5e5], # Particulate densities
    Kdet=1980.0, # Particulates detachment coefficient
    srcX=[(S, X, t, p) -> -b*X[1,:], # Source of particulates
          (S, X, t, p) ->  b*X[1,:]],
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> (mumax * S[1,:]) ./ (KM .+ S[1,:]) 
        (S, X, Lf, t, z, p) -> zeros(size(S[1,:]))],
    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Substrate"], # Substrate names
    Sin=[(t) -> 25],    # Substrate inflow (can be function of time)
    Sto=[25.0],          # Tank substrate concentration initial condition(s)
    Sbo=[0.0],          # Biofilm substrates concentration initial condition(s)
    Yxs=[0.378, 0],     # Biomass yield coefficient on substrate
    Daq=[1.38e-4],      # Substrate diffusion through boundary layer
    De=[6.9E-5],        # Substrate diffusion through biofilm     
    srcS=[(S,X,t,p) -> 0.0],     # Source of substrates
    
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
    Lfo=5.0e-6,     # Biofilm initial thickness [m]
    LL=1.0e-5,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,Xt,St,Pb,Sb,Lf,p) # Plot final results