using Biofilm 

# Define light as a function of time and depth within biofilm
diss=2000;  # Dissipation rate into biofilm [1/m]
smoothHeaviside(t,t0)=0.5*tanh.(100*(t.-t0).-0.5).+0.5
# Light :              turns off at t=0.25             turns on at t=0.75
intensity(t) = 1.0 - (smoothHeaviside(mod(t,1),0.25)-smoothHeaviside(mod(t,1),0.75))
# Dissipation of light into biofilm (1 at top with a rate of decrease of diss)
dissipation(z,Lf) = max.(0.0,1.0.-(Lf.-z)*diss)
light(t,z,Lf) = intensity(t)*dissipation(z,Lf)

# Source term constant
b = 0.1 

# Constants used for growthrates of particulate(s)
mumax = 0.4;    

# Tuple to hold all input parameters
p = (
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title =   "Phototroph Case",
    tFinal =  45,    # Simulation time [days]
    tol =     1e-4,  # Tolerance
    outPeriod = 5.0, # Time between outputs [days]
    # Let solver know when discontinuities (changes in light) occur
    discontinuityPeriod = 0.25,
    makePlots = false,

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames = ["Phototroph"],  # Particulate names
    Xto =  [1.0],             # Tank particulate concentration initial condition(s)
    Pbo =  [0.2],             # Biofilm particulates volume fraction initial condition(s) 
    rho =  [2.5e5],           # Particulate densities
    Kdet = 100.0,             # Particulates detachment coefficient
    srcX = [(S,X,Lf,t,z,p) -> 0.0], # Source of particulates
    # Growthrates for each particulate
    mu =[(S,X,Lf,t,z,p) -> mumax*light(t,z,Lf)],

    # ----------------- #
    # Solute Parameters #
    # ----------------- #
    SNames = ["Oxygen"],   # Solute names
    Sin =  [(t) -> 8.6],   # Solute inflow (can be function of time)
    Sto =  [8.6],          # Tank solute concentration initial condition(s)
    Sbo =  [8.6],          # Biofilm solutes concentration initial condition(s)
    Yxs =  [-0.52],        # Biomass yield coefficient on solute
    Dt =   [1.51e-4],      # Aquious solute diffusion through tank fluid
    Db =   [6.8E-5],       # Effective solute diffusion through biofilm
    srcS = [(S,X,Lf,t,z,p) -> 0.0],     # Source of solutes

    # --------------- #
    # Tank Parameters #
    # --------------- #
    V =0.01,       # Volume of tank [m³]
    A =  1,        # Surface area of biofilm [m²]
    Q = 10,        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz = 50,       # Number of grid points in biofilm
    Lfo =5.0e-6,   # Biofilm initial thickness [m]
    LL = 2.0e-4,  # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver