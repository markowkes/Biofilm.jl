using Biofilm 

# Constants used for growthrates of particulate(s)
mumax = 4;

# Source term constants
b=0.0

# Define light as a function of time and depth within biofilm
I=0.5;      # Light intesity -- Not needed (divides out)!!!
diss=1000;  # Dissipation rate into biofilm [1/m]
smoothHeaviside(t,t0)=tanh.(100*(t.-t0).-0.5)
# Light :         turns on at t=0.25             turns off at t=0.75
intensity(t) = I*(smoothHeaviside(mod(t,1),0.25)-smoothHeaviside(mod(t,1),0.75))
# Dissipation of light into biofilm (1 at top with a rate of decrease of diss)
dissipation(z,Lf) = max.(0.0,1.0.-(Lf.-z)*diss)
light(t,z,Lf) = intensity(t)*dissipation(z,Lf)


# Define a structure to hold all the parameters
p = param(
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> mumax*light(t,z,Lf)/I],

    # Source of particulates (constants defined above)
    src=[(S, X, p) -> 0.0],

    # Substrate inflow (can be function of time)
    Sin=[(t) -> 8.6],

    # Time
    tFinal=15,   # Simulation time [days]
    outPeriod=0.25,  # Time between outputs [days]

    # Simulation
    Title="Phototroph Case",
    SNames=["Oxygen"],
    XNames=["Phtotroph"],
    makePlots=true,

    # Tank Geometry
    V=0.01,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=1,          # Flowrate through tank [m³/s]
    Xo=[1.0],# Tank particulate initial condition(s)
    So=[8.6],    # Tank substrate initial condition(s)
    LL=2.0e-4,    # Boundary layer thickness [m]

    # Biofilm
    Nz=100,            # Number of grid points to represent biofilm
    Pbo=[0.2],     # Biofilm particulates initial condition(s)
    Sbo=[8.6],     # Biofilm substrates initial condition(s)
    Lfo=5.0E-6,    # Biofilm initial thickness [m]

    # Substance Constants
    Yxs=[0.0],     # Biomass yield coeffficient on substrate
    Daq=[1.51e-4],    # Substrate diffusion through boundary layer
    De =[6.8e-5],    # Substrate diffusion through biofilm     
    rho=[2.5e5],     # Particulate densities
    Kdet=100.0,     # Particulates detachment coefficient

    # Tolerance
    tol=1e-4,
)

sol,t,X,S,Pb,Sb,Lf = BiofilmSolver(p) # Run solver
outputs(t,X,S,Pb,Sb,Lf,p) # Plot final results