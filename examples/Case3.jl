using Biofilm 

# Create empty dictionary to hold parameters 
d = createDict()

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "Multiple Particulate Case")
addParam!(d, "tFinal",   100)   # Simulation time [days]
addParam!(d, "tol",      1e-4)  # Tolerance
addParam!(d, "outPeriod",5.0)   # Time between outputs [days]
addParam!(d, "optionalPlot","source") # 6th plot: "growthrate" (default) or "source"

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Living Bug","Dead Bug"])    # Particulate names
addParam!(d, "Xto",   [10.0, 0.0])        # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.08, 0.0])        # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.0e5, 2.0e5])     # Particulate densities
addParam!(d, "Kdet",  1980.0)             # Particulates detachment coefficient
b = 0.1 # Source term constant
addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> -b*X[1], # Source of particulates
                       (S,X,Lf,t,z,p) ->  b*X[1]])
# Growthrates for each particulate
mumax = 2; KM = 1;
addParam!(d, "mu", [(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1]) 
                    (S,X,Lf,t,z,p) -> 0.0 ] )

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Substrate"])     # Substrate names
addParam!(d, "Sin",   [(t) -> 25])   # Substrate inflow (can be function of time)
addParam!(d, "Sto",   [25.0])         # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [0.0])          # Biofilm substrates concentration initial condition(s)
addParam!(d, "Yxs",   [0.378, 0])        # Biomass yield coefficient on substrate
addParam!(d, "Daq",   [1.38e-4])       # Substrate diffusion through boundary layer
addParam!(d, "De",    [6.9E-5])       # Substrate diffusion through biofilm     
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
addParam!(d, "Lfo", 5.0e-6)   # Biofilm initial thickness [m]
addParam!(d, "LL",  1.0e-5)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
biofilm_plot(sol,p,size=(900,600))

# Save output
savefig("Case3.pdf")