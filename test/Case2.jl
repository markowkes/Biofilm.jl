using Biofilm

# Create empty dictionary to hold parameters 
d = createDict()

# Constants used for growthrates of particulate(s)
mumax = 20;
KM = 3;

# Define a structure to hold all the parameters
# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d,"Title","Multiple Independent Substrate Case")
addParam!(d,"tFinal",1)        # Simulation time [days]
addParam!(d,"tol",1e-2)        # Tolerance
addParam!(d,"outPeriod",0.1)   # Time between outputs [days]
addParam!(d,"makePlots",false)

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d,"XNames",["Bug"])     # Particulate names
addParam!(d,"Xto",[10.0])          # Tank particulate concentration initial condition(s)
addParam!(d,"Pbo",[0.2])          # Biofilm particulates volume fraction initial condition(s) 
addParam!(d,"rho",[1.0e5])        # Particulate densities
addParam!(d,"Kdet",1900.0)        # Particulates detachment coefficient
addParam!(d,"srcX",[(S,X,Lf,t,z,p) -> 0.0])     # Source of particulates
# Growthrates for each particulate (constants defined above!)
addParam!(d,"mu",[(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1])])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d,"SNames",["Substrate 1","Substrate 2"])  # Substrate names
addParam!(d,"Sin",[(t) -> 25, (t) -> 25]) # Substrate inflow (can be function of time)
addParam!(d,"Sto",[25.0,25.0])        # Tank substrate concentration initial condition(s)
addParam!(d,"Sbo",[0.0,0.0])          # Biofilm substrates concentration initial condition(s)
addParam!(d,"Yxs",[0.5 0.0])          # Biomass yield coefficient on substrate
addParam!(d,"Daq",[4.0e-5, 6.0e-5])   # Substrate diffusion through boundary layer
addParam!(d,"De",[1.0e-5, 1.5e-5])   # Substrate diffusion through biofilm     
addParam!(d,"srcS",[(S,X,Lf,t,z,p) -> 0.0,
                    (S,X,Lf,t,z,p) -> 0.0])     # Source of substrates
# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d,"V",0.1)        # Volume of tank [m³]
addParam!(d,"A",1)          # Surface area of biofilm [m²]
addParam!(d,"Q",1)          # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d,"Nz",50)          # Number of grid points in biofilm
addParam!(d,"Lfo",5.0e-6)     # Biofilm initial thickness [m]
addParam!(d,"LL",1.0e-4)      # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver