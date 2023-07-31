using Biofilm 

# Growthrate constants
KmB1 = 0.200; KmB3 = 11;  KmC2 = 20;  KI = 1.0;
mumaxA = 0.4;  mumaxB = 0.672;  mumaxC = 1.46;    

# Tuple to hold all input parameters
p = (
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title =   "SRB Test Case",
    tFinal =  500,   # Simulation time [days]
    tol =     1e-4,  # Tolerance
    outPeriod =10.0, # Time between outputs [days]
    plotPeriod = 20.0, # Time between outputs [days]

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames = ["SRB - Sulfate-Reducer"],    # Particulate names
    Xto =  [1.0e-6],  # Tank particulate concentration initial condition(s)
    Pbo =  [0.2],     # Biofilm particulates volume fraction initial condition(s) 
    rho =  [2.5e5],   # Particulate densities
    Kdet = 50.0,      # Particulates detachment coefficient
    srcX = [(S,X,Lf,t,z,p) -> 0.0], # Source of particulates
    # Growthrates for each particulate
    mu =[
        (S,X,Lf,t,z,p) ->  mumaxC*(S[2]./(KmC2.+S[2])).*(1.0./(1.0.+S[1]/KI))], #.-0.01*S[1]], # SRB

    # ----------------- #
    # Solute Parameters #
    # ----------------- #
    SNames = ["Oxygen","Sulfate","Hydrogen Sulfide"],     # Solute names
    Sin =  [
        (t) -> 8.6    # Solute inflow (can be function of time)
        (t) -> 48.0
        (t) -> 0.0 ],
    Sto =  [8.6,48.0,0.0],          # Tank substrate concentration initial condition(s)
    Sbo =  [8.6,48.0,1e-5],         # Biofilm substrates concentration initial condition(s)
    # Biomass yield coefficient on substrate
    Yxs =  [ # oxygen  sulfate  Hy. sulfide
                0.00    0.584   -1.645],    # SRB uses sulfate and produces sulfide
    Dt =   [1.51e-4,8e-5,1.21e-4],   # Aquious substrate diffusion through tank fluid
    Db =   [6.8e-5,4e-5,6.04e-5],    # Effective substrate diffusion through biofilm
    srcS = [
        (S,X,Lf,t,z,p) -> 0.0,       # Source of substrates
        (S,X,Lf,t,z,p) -> 0.0,
        (S,X,Lf,t,z,p) -> 0.0],

    # --------------- #
    # Tank Parameters #
    # --------------- #
    V =0.1,        # Volume of tank [m³]
    A =  1,        # Surface area of biofilm [m²]
    Q = 10,        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz = 20,       # Number of grid points in biofilm
    Lfo =5.0e-6,   # Biofilm initial thickness [m]
    LL = 2.0e-4,  # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver