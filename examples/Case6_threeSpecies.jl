using Biofilm 

# Constants used for growthrates of particulate(s)
KmB1 = 0.200; KmB3 = 11;  KmC2 = 20;  KI = 0.5;
mumaxA = 0.4;  mumaxB = 0.672;  mumaxC = 1.46;

## Define light as a function of time and depth within biofilm
smoothHeaviside(t,t0)=0.5*tanh.(100*(t.-t0).-0.5).+0.5

# Light :              turns off at t=0.25             turns on at t=0.75
intensity(t) = 1.0 - (smoothHeaviside(mod(t,1),0.25)-smoothHeaviside(mod(t,1),0.75))

# Dissipation of light into biofilm (1 at top with a rate of decrease of diss)
diss=1500;  # Dissipation rate into biofilm [1/m]
dissipation(z,Lf) = max.(0.0,1.0.-(Lf.-z)*diss)

# Light at a given location 
light(t,z,Lf) = intensity(t)*dissipation(z,Lf)

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Ramsing et al. 1993 Test Case",
    tFinal=500,     # Simulation time [days]
    tol=1e-4,       # Tolerance
    outPeriod=5,    # Time between outputs [days]
    plotPeriod=10,    # Time between plots [days] (make multiple of outPeriod!)

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["Phototroph","SOB - Sulfide-Oxidizer","SRB - Sulfate-Reducer"], # Particulate names
    Xo=[1.0,1.0,1.0], # Tank particulate concentration initial condition(s)
    Pbo=[0.2/3,0.2/3,0.2/3], # Biofilm particulates volume fraction initial condition(s) 
    rho=[2.5e5,2.5e5,2.5e5], # Particulate densities
    Kdet=50.0, # Particulates detachment coefficient
    src=[(S, X, t, p) -> 0.0, # Source of particulates
         (S, X, t, p) -> 0.0,
         (S, X, t, p) -> 0.0],
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> mumaxA*light(t,z,Lf),
        (S, X, Lf, t, z, p) -> mumaxB*(S[1,:]./(KmB1.+S[1,:])).*(S[3,:]./(KmB3.+S[3,:])),
        (S, X, Lf, t, z, p) -> mumaxC*(S[2,:]./(KmC2.+S[2,:])).*(1.0./(1.0.+S[1,:]/KI))],
    discontinuityPeriod=0.25,  # Let solver know when discontinuities (changes in light) occur
    
    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Oxygen","Sulfate","Hydrogen Sulfide"], # Substrate names
    Sin=[(t) -> 8.6    # Substrate inflow (can be function of time)
         (t) -> 48.0
         (t) -> 0.0],
    So=[8.6,48.0,0.0],  # Tank substrate concentration initial condition(s)
    Sbo=[8.6,48.0,0.0], # Biofilm substrates concentration initial condition(s)
    # Biomass yield coefficient on substrate
    #     oxygen  sulfate  Hy. sulfide
    Yxs=[ -0.52    0.0      0.0        # Phototropho produces oxygen
           0.058   0.0      0.09       # SOB uses oxygen and sulfide
           0.00    0.584   -1.645],    # SRB uses sulfate and produces sulfide
    Daq=[1.51e-4,8e-5,1.21e-4],    # Substrate diffusion through boundary layer
    De =[6.8e-5,4e-5,6.04e-5],     # Substrate diffusion through biofilm     
           
    
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.01,       # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=10,         # Flowrate through tank [m³/s]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=20,          # Number of grid points in biofilm
    Lfo=5.0e-6,     # Biofilm initial thickness [m]
    LL=2.0e-4,      # Boundary layer thickness [m]
)

t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,X,S,Pb,Sb,Lf,p) # Plot final results