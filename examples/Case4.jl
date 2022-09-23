using Biofilm 

# Constants used for growthrates of particulate(s)
mumax = 2000;
KM = 2500;

# Source term constants
b=0.0

# Define a structure to hold all the parameters
p = param(
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> mumax * S[1,:] ./ KM 
        (S, X, Lf, t, z, p) -> mumax * S[2,:] ./ KM ],

    # Source of particulates (constants defined above)
    src=[(S, X, t, p) -> -b*X[1,:],
         (S, X, t, p) ->  b*X[1,:]],

    # Substrate inflow (can be function of time)
    Sin=[(t) -> 25,
         (t) -> 25],

    # Time
    tFinal=5,   # Simulation time [days]
    outPeriod=0.5,  # Time between outputs [days]

    # Simulation
    Title="Multiple Particulates and Substrates Case",
    SNames=["Substrate 1", "Substrate 2"],
    XNames=["Bug 1","Bug 2"],
    makePlots=true,

    # Tank Geometry
    V=0.1,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=1,          # Flowrate through tank [m³/s]
    Xo=[10.0,10.0],# Tank particulate initial condition(s)
    So=[25.0,25.0],    # Tank substrate initial condition(s)
    LL=1.0e-4,    # Boundary layer thickness [m]

    # Biofilm
    Nz=50,            # Number of grid points to represent biofilm
    Pbo=[0.2,0.2],     # Biofilm particulates initial condition(s)
    Sbo=[25.0,25.0],     # Biofilm substrates initial condition(s)
    Lfo=5.0E-6,    # Biofilm initial thickness [m]

    # Substance Constants
    Yxs=[0.5 0.0
         0.0 0.278],     # Biomass yield coeffficient on substrate
    Daq=[4.0e-5,6.0e-5],    # Substrate diffusion through boundary layer
    De =[1.0e-5,1.5e-5],    # Substrate diffusion through biofilm     
    rho=[3.0e5,3.0e5],     # Particulate densities
    Kdet=1900.0,     # Particulates detachment coefficient

    # Tolerance
    tol=1e-4,
)

t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,zm,X,S,Pb,Sb,Lf,p) # Plot final results