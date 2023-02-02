using Biofilm 

# Constants used for growthrates of particulate(s)
mumax = 0.4;

# Define light as a function of time and depth within biofilm
diss=2000;  # Dissipation rate into biofilm [1/m]
smoothHeaviside(t,t0)=0.5*tanh.(100*(t.-t0).-0.5).+0.5
# Light :              turns off at t=0.25             turns on at t=0.75
intensity(t) = 1.0 - (smoothHeaviside(mod(t,1),0.25)-smoothHeaviside(mod(t,1),0.75))
# Dissipation of light into biofilm (1 at top with a rate of decrease of diss)
dissipation(z,Lf) = max.(0.0,1.0.-(Lf.-z)*diss)
light(t,z,Lf) = intensity(t)*dissipation(z,Lf)


# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Phototroph Case",
    tFinal=45,     # Simulation time [days]
    tol=1e-4,       # Tolerance
    outPeriod=5,    # Time between outputs [days]

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["Phototroph"], # Particulate names
    Xto=[1.0],  # Tank particulate concentration initial condition(s)
    Pbo=[0.2], # Biofilm particulates volume fraction initial condition(s) 
    rho=[2.5e5], # Particulate densities
    Kdet=100.0, # Particulates detachment coefficient
    srcX=[(St, Xt, t, p) -> 0.0], # Source of particulates
    # Growthrates for each particulate (constants defined above!)
    mu=[(St, Xt, Lf, t, z, p) -> mumax*light(t,z,Lf)],
    discontinuityPeriod=0.25,  # Let solver know when discontinuities (changes in light) occur

    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Oxygen"], # Substrate names
    Sin=[(t) -> 8.6],    # Substrate inflow (can be function of time)
    Sto=[8.6],          # Tank substrate concentration initial condition(s)
    Sbo=[8.6],          # Biofilm substrates concentration initial condition(s)
    Yxs=[-0.52],     # Biomass yield coefficient on substrate
    Daq=[1.51e-4],      # Substrate diffusion through boundary layer
    De=[6.8E-5],        # Substrate diffusion through biofilm     
    srcS=[(St,Xt,t,p) -> 0.0],     # Source of substrates
    
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.01,       # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=10,         # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=50,          # Number of grid points in biofilm
    Lfo=5.0e-6,     # Biofilm initial thickness [m]
    LL=2.0e-4,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,Xt,St,Pb,Sb,Lf,p) # Plot final results