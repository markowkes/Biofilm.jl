using Biofilm 
using Plots

# Create empty dictionary to hold parameters 
d = createDict()

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "Single Substrate and Particulate Case")
addParam!(d, "tFinal",   1.0)   # Simulation time [days]
addParam!(d, "tol",      1e-2)  # Tolerance
addParam!(d, "outPeriod",0.1)   # Time between outputs [days]

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Aerobe"])    # Particulate names
addParam!(d, "Xto",   [10.0])        # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.08])        # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.0E4])       # Particulate densities
addParam!(d, "Kdet",  20000.0)       # Particulates detachment coefficient
addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> 0.0])     # Source of particulates
# Growthrates for each particulate
mumax = 20; KM = 3;
addParam!(d, "mu", [(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1])])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Oxygen"])     # Substrate names
addParam!(d, "Sin",   [(t) -> 100])   # Substrate inflow (can be function of time)
addParam!(d, "Sto",   [10.0])         # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [0.0])          # Biofilm substrates concentration initial condition(s)
addParam!(d, "Yxs",   [2.646])        # Biomass yield coefficient on substrate
addParam!(d, "Dt",    [4.0E-5])       # Aquious substrate diffusion through tank fluid
addParam!(d, "Db",    [6.9E-5])       # Effective substrate diffusion through biofilm
addParam!(d, "srcS",  [(S,X,Lf,t,z,p) -> 0.0])     # Source of substrates

# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d, "V", 0.1)        # Volume of tank [m³]
addParam!(d, "A",   1)        # Surface area of biofilm [m²]
addParam!(d, "Q",   1)        # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d, "Nz",  50)       # Number of grid points in biofilm
addParam!(d, "Lfo", 1.0E-5)   # Biofilm initial thickness [m]
addParam!(d, "LL",  1.00E-7)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
biofilm_plot(sol,p,size=(900,600))

# Save output
savefig("Case1.pdf")
