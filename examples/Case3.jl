using Biofilm 

# Constants used for growthrates of particulate(s)
mumax = 2;
KM = 1;

# Source term constants
b=0.1

# Define a structure to hold all the parameters
p = param(
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> (mumax * S[1,:]) ./ (KM .+ S[1,:]) 
        (S, X, Lf, t, z, p) -> zeros(size(S[1,:]))],

    # Source of particulates (constants defined above)
    src=[(S, X, p) -> -b*X[1,:],
         (S, X, p) ->  b*X[1,:]],

    # Substrate inflow (can be function of time)
    Sin=[(t) -> 25],

    # Time
    tFinal=100,   # Simulation time [days]
    outPeriod=5,  # Time between outputs [days]

    # Simulation
    Title="Multiple Particulate Case",
    SNames=["Substrate"],
    XNames=["Bug 1","Bug 2"],
    makePlots=true,

    # Tank Geometry
    V=0.1,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=1,          # Flowrate through tank [m³/s]
    Xo=[10.0,0.0],# Tank particulate initial condition(s)
    So=[25.0],    # Tank substrate initial condition(s)
    LL=1.0e-5,    # Boundary layer thickness [m]

    # Biofilm
    Nz=50,            # Number of grid points to represent biofilm
    Pbo=[0.08,0.0],     # Biofilm particulates initial condition(s)
    Sbo=[0.0],     # Biofilm substrates initial condition(s)
    Lfo=5.0E-6,    # Biofilm initial thickness [m]

    # Substance Constants
    Yxs=[0.378, 0],     # Biomass yield coeffficient on substrate
    Daq=[1.38e-4],    # Substrate diffusion through boundary layer
    De =[6.9e-5],    # Substrate diffusion through biofilm     
    rho=[2.5e5,2.5e5],     # Particulate densities
    Kdet=1980.0,     # Particulates detachment coefficient

    # Tolerance
    tol=1e-4,
)

sol,t,X,S,Pb,Sb,Lf = BiofilmSolver(p) # Run solver
outputs(t,X,S,Pb,Sb,Lf,p) # Plot final results