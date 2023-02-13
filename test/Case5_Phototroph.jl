using Biofilm 

# Create empty dictionary to hold parameters 
d = createDict()

# Define light as a function of time and depth within biofilm
diss=2000;  # Dissipation rate into biofilm [1/m]
smoothHeaviside(t,t0)=0.5*tanh.(100*(t.-t0).-0.5).+0.5
# Light :              turns off at t=0.25             turns on at t=0.75
intensity(t) = 1.0 - (smoothHeaviside(mod(t,1),0.25)-smoothHeaviside(mod(t,1),0.75))
# Dissipation of light into biofilm (1 at top with a rate of decrease of diss)
dissipation(z,Lf) = max.(0.0,1.0.-(Lf.-z)*diss)
light(t,z,Lf) = intensity(t)*dissipation(z,Lf)

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "Phototroph Case")
addParam!(d, "tFinal",   45)   # Simulation time [days]
addParam!(d, "tol",      1e-4)  # Tolerance
addParam!(d, "outPeriod",5.0)   # Time between outputs [days]
# Let solver know when discontinuities (changes in light) occur
addParam!(d, "discontinuityPeriod",0.25)  
addParam!(d, "makePlots",false)

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Phototroph"])    # Particulate names
addParam!(d, "Xto",   [1.0])             # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.2])             # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.5e5])           # Particulate densities
addParam!(d, "Kdet",  100.0)             # Particulates detachment coefficient
b = 0.1 # Source term constant
addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> 0.0]) # Source of particulates
# Growthrates for each particulate
# Constants used for growthrates of particulate(s)
mumax = 0.4;
addParam!(d, "mu", [(S,X,Lf,t,z,p) -> mumax*light(t,z,Lf)])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Oxygen"])     # Substrate names
addParam!(d, "Sin",   [(t) -> 8.6])   # Substrate inflow (can be function of time)
addParam!(d, "Sto",   [8.6])          # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [8.6])          # Biofilm substrates concentration initial condition(s)
addParam!(d, "Yxs",   [-0.52])        # Biomass yield coefficient on substrate
addParam!(d, "Daq",   [1.51e-4])       # Substrate diffusion through boundary layer
addParam!(d, "De",    [6.8E-5])       # Substrate diffusion through biofilm     
addParam!(d, "srcS",  [(S,X,Lf,t,z,p) -> 0.0])     # Source of substrates

# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d, "V", 0.01)       # Volume of tank [m³]
addParam!(d, "A",   1)        # Surface area of biofilm [m²]
addParam!(d, "Q",  10)        # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d, "Nz",  50)       # Number of grid points in biofilm
addParam!(d, "Lfo", 5.0e-6)   # Biofilm initial thickness [m]
addParam!(d, "LL",  2.0e-4)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver