using Biofilm 

# Constants used for growthrates of particulate(s)
KmB1 = 0.200; KmB3 = 11;  KmC2 = 20;  KI = 1.0;
mumaxA = 0.4;  mumaxB = 0.672;  mumaxC = 1.46;

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="SOB-SRB Test Case",
    tFinal=500,     # Simulation time [days]
    tol=1e-4,       # Tolerance
    outPeriod=10,   # Time between outputs [days]
    plotPeriod=20,  # Time between plots [days] (make multiple of outPeriod!)
    #optionalPlot="source", # 6th plot: "growthrate" (default) or "source"
    
    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["SOB - Sulfide-Oxidizer","SRB - Sulfate-Reducer"], # Particulate names
    Xto=[1.0e-6,1.0e-6],  # Tank particulate concentration initial condition(s)
    Pbo=[0.2/2,0.2/2], # Biofilm particulates volume fraction initial condition(s) 
    rho=[2.5e5,2.5e5], # Particulate densities
    Kdet=50.0, # Particulates detachment coefficient
    srcX=[(St, Xt, t, p) -> 0.0, # Source of particulates
         (St, Xt, t, p) -> -10*S[1,:]], # dies near oxygen
    # Growthrates for each particulate (constants defined above!)
    mu=[(St, Xt, Lf, t, z, p) -> mumaxB*(S[1,:]./(KmB1.+S[1,:])).*(S[3,:]./(KmB3.+S[3,:])), # SOB
        (St, Xt, Lf, t, z, p) -> mumaxC*(S[2,:]./(KmC2.+S[2,:])).*(1.0./(1.0.+S[1,:]/KI))], # SRB

    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Oxygen","Sulfate","Hydrogen Sulfide"], # Substrate names
    Sin=[(t) -> 8.6     # Substrate inflow (can be function of time)
         (t) -> 48.0
         (t) -> 0.0],
    Sto=[8.6,48.0,0.0],   # Tank substrate concentration initial condition(s)
    Sbo=[8.6,48.0,1e-5], # Biofilm substrates concentration initial condition(s)
    # Biomass yield coefficient on substrate
    #     oxygen  sulfate  Hy. sulfide
    Yxs=[  0.058    0.0      0.09       # SOB uses oxygen and sulfide
           0.00    0.584   -1.645],    # SRB uses sulfate and produces sulfide
    Daq=[1.51e-4,8e-5,1.21e-4],    # Substrate diffusion through boundary layer
    De =[6.8e-5,4e-5,6.04e-5],     # Substrate diffusion through biofilm     
    srcS=[(St,Xt,t,p) -> 0.0,
          (St,Xt,t,p) -> 0.0,       # Source of substrates
          (St,Xt,t,p) -> 0.0],
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.01,       # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=10,         # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=20,          # Number of grid points in biofilm
    Lfo=5.0e-6,     # Biofilm initial thickness [m]
    LL=2.0e-4,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,Xt,St,Pb,Sb,Lf,p) # Plot final results