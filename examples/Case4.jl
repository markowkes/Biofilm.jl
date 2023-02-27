using Biofilm 

# Create empty dictionary to hold parameters 
d = createDict()

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "SOB-SRB Test Case")
addParam!(d, "tFinal",   100)    # Simulation time [days]
addParam!(d, "tol",      1e-4)   # Tolerance
addParam!(d, "outPeriod", 1.0)   # Time between outputs [days]
addParam!(d, "plotPeriod",5.0)   # Time between outputs [days]

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["SOB - Sulfide-Oxidizer","SRB - Sulfate-Reducer","Dead Bacteria"])    # Particulate names
addParam!(d, "Xto",   [1.0e-6,1.0e-6,0.0])      # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.2/2,0.2/2,0.0])        # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.5e5,2.5e5,2.5e5])      # Particulate densities
addParam!(d, "Kdet",  500.0)                    # Particulates detachment coefficient
# Source for each particulate
D_O_SRB = 10.0
D_SOB = 1e-2
D_SRB = 1e-1
addParam!(d, "srcX",  [
    (S,X,Lf,t,z,p) -> - D_SOB*X[1]                             ,  # SOB - slowly dies
    (S,X,Lf,t,z,p) ->              - D_SRB*X[2] - D_O_SRB*S[1] ,  # SRB dies near oxygen and slowly dies
    (S,X,Lf,t,z,p) -> + D_SRB*X[2] + D_SOB*X[1] + D_O_SRB*S[1] ]) # Dead bacteria (opposite of above) 
# Growthrates for each particulate
KmB1 = 0.200; KmB3 = 11;  KmC2 = 20;  KI = 1.0;
mumaxA = 0.4;  mumaxB = 0.672;  mumaxC = 10.46;
addParam!(d, "mu", [
    (S,X,Lf,t,z,p) -> mumaxB*(S[1]./(KmB1.+S[1])).*(S[3]./(KmB3.+S[3])),   # SOB
    (S,X,Lf,t,z,p) -> mumaxC*(S[2]./(KmC2.+S[2])).*(1.0./(1.0.+S[1]/KI)) , # SRB
    (S,X,Lf,t,z,p) -> 0.0 ] )                                              # Dead

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Oxygen","Sulfate","Hydrogen Sulfide"])     # Substrate names
addParam!(d, "Sin",   [
    (t) -> 8.6     # Substrate inflow (can be function of time)
    (t) -> 48.0
    (t) -> 0.0] ) 
addParam!(d, "Sto",   [8.6,48.0,0.0])         # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [8.6,48.0,1e-5])          # Biofilm substrates concentration initial condition(s)
# Biomass yield coefficient on substrate
    
addParam!(d, "Yxs",   [#oxygen  sulfate  Hy. sulfide
                        0.058    0.0      0.09    # SOB uses oxygen and sulfide
                        0.00    0.584   -1.645    # SRB uses sulfate and produces sulfide
                        0.00     0.0      0.0 ])  # Dead not needed/used for growth
addParam!(d, "Daq",   [1.51e-4,8e-5,1.21e-4])     # Substrate diffusion through boundary layer
addParam!(d, "De",    [6.8e-5,4e-5,6.04e-5])      # Substrate diffusion through biofilm     
addParam!(d, "srcS",  [(S,X,Lf,t,z,p) -> 0.0,     # Source of substrates
                       (S,X,Lf,t,z,p) -> 0.0,
                       (S,X,Lf,t,z,p) -> 0.0])     

# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d, "V", 0.1)        # Volume of tank [m³]
addParam!(d, "A",   1)        # Surface area of biofilm [m²]
addParam!(d, "Q",  10)        # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d, "Nz",  40)       # Number of grid points in biofilm
addParam!(d, "Lfo", 5.0e-6)   # Biofilm initial thickness [m]
addParam!(d, "LL",  2.0e-4)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
biofilm_plot(sol,p)