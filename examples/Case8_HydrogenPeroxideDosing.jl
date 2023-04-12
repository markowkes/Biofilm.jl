using Biofilm 

# Create empty dictionary to hold parameters 
d = createDict()

smoothHeaviside(t,t0)=0.5*tanh.(10*(t.-t0).-0.5).+0.5

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "Hydrogen Peroxide Dosing")
addParam!(d, "tFinal",   10)    # Simulation time [days]
addParam!(d, "tol",      1e-8)  # Tolerance
addParam!(d, "outPeriod",1.0)   # Time between outputs [days]
addParam!(d, "discontinuityPeriod",2.5) # Let solver know when discontinuities

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Live","Dead"])    # Particulate names
addParam!(d, "Xto",   [1.0, 0.0])         # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.08, 0.0])        # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.5e5, 2.5e5])     # Particulate densities
addParam!(d, "Kdet",  1e4)                # Particulates detachment coefficient
k_dis = 0.5; #m³/g/d   # Source term constant
addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> -k_dis*X[1].*S[2], # Source of particulates
                       (S,X,Lf,t,z,p) -> +k_dis*X[1].*S[2] ])
# Growthrates for each particulate
mumax = 9.6; # 1/days
KM = 5; # g/m³
addParam!(d, "mu", [(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1]) 
                    (S,X,Lf,t,z,p) -> 0.0 ] )

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Glucose","Hydrogen Peroxide"])     # Substrate names
addParam!(d, "Sin",   [
        (t) -> 100,    # Substrate inflow (can be function of time)
        (t) -> 500*smoothHeaviside(t,2.5)])
addParam!(d, "Sto",   [100.0, 0.0])          # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [  0.0, 0.0])          # Biofilm substrates concentration initial condition(s)
# Biomass yield coefficient on substrate
addParam!(d, "Yxs",   [#Glucose   H. Per.
                        0.26       0.0       # Live use glucose
                        0.00       0.0   ])  # Dead doesn't use/produce anything
addParam!(d, "Dt",    [5.2e-5, 1.09e-4])     # Aquious substrate diffusion through tank fluid
addParam!(d, "Db",    [1.3e-5, 6.52e-5])     # Effective substrate diffusion through biofilm
k_b  = 10.0; #m³/g/d
addParam!(d, "srcS",  [                      # Source of substrates
    (S,X,Lf,t,z,p) -> 0.0,  
    (S,X,Lf,t,z,p) -> -k_b*X[1].*S[2].-k_b*X[2].*S[2] ])

# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d, "V", 0.1)        # Volume of tank [m³]
addParam!(d, "A",   1)        # Surface area of biofilm [m²]
addParam!(d, "Q",   1)        # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d, "Nz",  50)        # Number of grid points in biofilm
addParam!(d, "Lfo", 50.0e-6)   # Biofilm initial thickness [m]
addParam!(d, "LL",  1.0e-5)    # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver