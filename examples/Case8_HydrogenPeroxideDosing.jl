using Biofilm 

# Constants used for growthrates of particulate(s)
Km = 5; # g/m³
mumax = 9.6; # 1/days
k_dis = 0.5; #m³/g/d
k_b  = 10.0; #m³/g/d

smoothHeaviside(t,t0)=0.5*tanh.(10*(t.-t0).-0.5).+0.5

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Hydrogen Peroxide Dosing",
    tFinal=10,     # Simulation time [days]
    tol=1e-8,       # Tolerance
    outPeriod=1,    # Time between outputs [days]
    plotPeriod=1,    # Time between plots [days] (make multiple of outPeriod!)
    #optionalPlot="source", # 6th plot: "growthrate" (default) or "source"
    discontinuityPeriod=2.5,  # Let solver know when discontinuities (changes in light) occur

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["Live","Dead"], # Particulate names
    Xto=[1.0,0.0], # Tank particulate concentration initial condition(s)
    Pbo=[0.08,0.0], # Biofilm particulates volume fraction initial condition(s) 
    rho=[2.5e5,2.5e5], # Particulate densities
    Kdet=1e4, # Particulates detachment coefficient
    srcX=[(S, X, t, p) -> -k_dis*X[1,:].*S[2,:], 
          (S, X, t, p) -> +k_dis*X[1,:].*S[2,:]],
    # Growthrates for each particulate (constants defined above!)
    mu=[(S,X,Lf,t,z,p) -> mumax * S[1,:] ./ (Km .+ S[1,:]), 
        (S,X,Lf,t,z,p) -> zeros(size(S[1,:])) ],
    
    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Glucose","Hydrogen Peroxide"], # Substrate names
    Sin=[(t) -> 100,    # Substrate inflow (can be function of time)
         (t) -> 500*smoothHeaviside(t,2.5)],
    Sto=[100.0,0.0],  # Tank substrate concentration initial condition(s)
    Sbo=[0.0,0.0], # Biofilm substrates concentration initial condition(s)
    # Biomass yield coefficient on substrate
    #     Glucose   H. Per.
    Yxs=[ 0.26       0.0        # Live use glucose
          0.00       0.0   ],   # Dead doesn't use/produce anything
    Daq=[5.2e-5, 1.09e-4],    # Substrate diffusion through boundary layer
    De =[1.3e-5, 6.52e-5],     # Substrate diffusion through biofilm     
    srcS=[(S, X, t, p) -> 0.0,  
          (S, X, t, p) -> -k_b*X[1,:].*S[2,:].-k_b*X[2,:].*S[2,:] ],
           
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.1,        # Volume of tank [m³]
    A=1.0,        # Surface area of biofilm [m²]
    Q=1.0,        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=50,          # Number of grid points in biofilm
    Lfo=50.0e-6,     # Biofilm initial thickness [m]
    LL=1.0e-5,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,Xt,St,Pb,Sb,Lf,p) # Plot final results