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
addParam!(d, "Title",    "Ramsing et al. 1993 Test Case")
addParam!(d, "tFinal",   500)   # Simulation time [days]
addParam!(d, "tol",      1e-4)  # Tolerance
addParam!(d, "outPeriod",  5.0) # Time between outputs [days]
addParam!(d, "plotPeriod",10.0) # Time between outputs [days]
# Let solver know when discontinuities (changes in light) occur
addParam!(d, "discontinuityPeriod",0.25)  


# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Phototroph","SOB - Sulfide-Oxidizer","SRB - Sulfate-Reducer"])    # Particulate names
addParam!(d, "Xto",   [1.0,1.0,1.0])        # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.2/3,0.2/3,0.2/3])  # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.5e5,2.5e5,2.5e5])  # Particulate densities
addParam!(d, "Kdet",  50.0)                 # Particulates detachment coefficient
addParam!(d, "srcX",  [
    (S,X,Lf,t,z,p) -> 0.0, # Source of particulates
    (S,X,Lf,t,z,p) -> 0.0,
    (S,X,Lf,t,z,p) -> 0.0])
# Growthrates for each particulate
KmB1 = 0.200; KmB3 = 11;  KmC2 = 20;  KI = 0.5;
mumaxA = 0.4;  mumaxB = 0.672;  mumaxC = 1.46;
addParam!(d, "mu", [
    (S,X,Lf,t,z,p) -> mumaxA*light(t,z,Lf),
    (S,X,Lf,t,z,p) -> mumaxB*(S[1]./(KmB1.+S[1])).*(S[3]./(KmB3.+S[3])),
    (S,X,Lf,t,z,p) -> mumaxC*(S[2]./(KmC2.+S[2])).*(1.0./(1.0.+S[1]/KI))])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Oxygen","Sulfate","Hydrogen Sulfide"])     # Substrate names
addParam!(d, "Sin",   [
    (t) -> 8.6    # Substrate inflow (can be function of time)
    (t) -> 48.0
    (t) -> 0.0 ])
addParam!(d, "Sto",   [8.6,48.0,0.0])          # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [8.6,48.0,0.0])          # Biofilm substrates concentration initial condition(s)
# Biomass yield coefficient on substrate
addParam!(d, "Yxs",   [ # oxygen  sulfate  Hy. sulfide
                          -0.52    0.0      0.0       # Phototropho produces oxygen
                          0.058   0.0      0.09       # SOB uses oxygen and sulfide
                          0.00    0.584   -1.645])    # SRB uses sulfate and produces sulfide
addParam!(d, "Daq",   [1.51e-4,8e-5,1.21e-4])       # Substrate diffusion through boundary layer
addParam!(d, "De",    [6.8e-5,4e-5,6.04e-5])       # Substrate diffusion through biofilm     
addParam!(d, "srcS",  [
    (S,X,Lf,t,z,p) -> 0.0,       # Source of substrates
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
addParam!(d, "Nz",  20)       # Number of grid points in biofilm
addParam!(d, "Lfo", 5.0e-6)   # Biofilm initial thickness [m]
addParam!(d, "LL",  2.0e-4)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver