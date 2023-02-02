using Biofilm 

# Constants used for growthrates of particulate(s)
mumax = 2000;
KM = 2500;

# Source term constants
b=0.0

# Define a structure to hold all the parameters
p = param(
     # --------------------- #
     # Simulation Parameters #
     # --------------------- #
     Title="Multiple Particulates and Substrates Case",
     tFinal=5,      # Simulation time [days]
     tol=1e-4,      # Tolerance
     outPeriod=0.5, # Time between outputs [days]
     makePlots=false,

     # ---------------------- #
     # Particulate Parameters #
     # ---------------------- #
     XNames=["Bug 1","Bug 2"], # Particulate names
     Xto=[10.0,10.0],  # Tank particulate concentration initial condition(s)
     Pbo=[0.2,0.2], # Biofilm particulates volume fraction initial condition(s) 
     rho=[3e5,3e5], # Particulate densities
     Kdet=1900.0, # Particulates detachment coefficient
     srcX=[(St, Xt, t, p) -> -b*X[1,:], # Source of particulates
          (St, Xt, t, p) ->  b*X[1,:]],
     # Growthrates for each particulate (constants defined above!)
     mu=[(St, Xt, Lf, t, z, p) -> mumax * S[1,:] ./ KM 
         (St, Xt, Lf, t, z, p) -> mumax * S[2,:] ./ KM ],

     # -------------------- #
     # Substrate Parameters #
     # -------------------- #
     SNames=["Substrate 1", "Substrate 2"], # Substrate names
     Sin=[(t) -> 25,     # Substrate inflow (can be function of time)
          (t) -> 25],
     Sto=[25.0,25.0],     # Tank substrate concentration initial condition(s)
     Sbo=[25.0,25.0],    # Biofilm substrates concentration initial condition(s)
     Yxs=[0.5  0.0       # Biomass yield coefficient on substrate
          0.0  0.278],  
     Daq=[4.0e-5,6.0e-5], # Substrate diffusion through boundary layer
     De=[1.0e-5,1.5e-5],  # Substrate diffusion through biofilm     
     srcS=[(St,Xt,t,p) -> 0.0,
          (St,Xt,t,p) -> 0.0],     # Source of substrates
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
     LL=1.0e-4,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
