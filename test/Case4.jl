using Biofilm 

# Create empty dictionary to hold parameters 
d = createDict()

# Constants used for growthrates of particulate(s)
mumax = 2000;
KM = 2500;

# Source term constants
b=0.0

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d,"Title","Multiple Particulates and Substrates Case")
addParam!(d,"tFinal",5)      # Simulation time [days]
addParam!(d,"tol",1e-4)      # Tolerance
addParam!(d,"outPeriod",0.5) # Time between outputs [days]
addParam!(d,"makePlots",false)

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d,"XNames",["Bug 1","Bug 2"]) # Particulate names
addParam!(d,"Xto",[10.0,10.0])  # Tank particulate concentration initial condition(s)
addParam!(d,"Pbo",[0.2,0.2]) # Biofilm particulates volume fraction initial condition(s) 
addParam!(d,"rho",[3e5,3e5]) # Particulate densities
addParam!(d,"Kdet",1900.0) # Particulates detachment coefficient
addParam!(d,"srcX",[(S,X,Lf,t,z,p) -> -b*X[1]  # Source of particulates
                    (S,X,Lf,t,z,p) ->  b*X[1] ])
# Growthrates for each particulate (constants defined above!)
addParam!(d,"mu",[(S,X,Lf,t,z,p) -> mumax * S[1] ./ KM 
                  (S,X,Lf,t,z,p) -> mumax * S[2] ./ KM ])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d,"SNames",["Substrate 1", "Substrate 2"]) # Substrate names
addParam!(d,"Sin",[(t) -> 25,     # Substrate inflow (can be function of time)
                   (t) -> 25])
addParam!(d,"Sto",[25.0,25.0])     # Tank substrate concentration initial condition(s)
addParam!(d,"Sbo",[25.0,25.0])    # Biofilm substrates concentration initial condition(s)
addParam!(d,"Yxs",[0.5  0.0       # Biomass yield coefficient on substrate
                   0.0  0.278])
addParam!(d,"Daq",[4.0e-5,6.0e-5]) # Substrate diffusion through boundary layer
addParam!(d,"De",[1.0e-5,1.5e-5])  # Substrate diffusion through biofilm     
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