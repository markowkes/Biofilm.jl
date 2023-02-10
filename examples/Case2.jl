using Biofilm 
using Plots

# This case is adapted from the Acid Stress Response simulation in https://doi.org/10.1371/journal.pone.0083626

# Notation
# s = Glucose
# p = Lactate

# Constants from paper (with unit conversions)
Ds  = 2.24e-10 * 86400 # Glucose diffusion coefficient [m²/s -> m²/d]
Dp  = 2.24e-10 * 86400 # Lactate diffusion coefficient [m²/s -> m²/d]
L₀  = 80.0 * 1e-6 # Initial biofilm thickness [μm -> mm]
Rp  = 100.0 * 1e-6 # Penetration depth of glucose [μm -> mm]
X   = 5.0 * 1000 # Cell density [g/m³]
Yxs = 0.5 # Yield coefficient of biomass on glucose g_X/g_s
Yps = 0.9 # Yield coefficient of lactate on glucose g_p/g_s 
h₀  = 100 * 1e-6 # Characteristic length scale [μm -> m]
m   = 0.1 * 24 # Turnover rate of gene [1/h -> 1/d]
s⁰  = 800  # Bulk concentration of glucose [mg/l = g/m³]
t₀  = 3600 / 86400 # Characteristic time scale [s -> d]
μᵢ  = 0.001 * 24 # Specific growth rate coefficient [1/h⋅l/mg -> 1/d⋅m³/g]
μₚ  = 3000.0 # Cutoff used in growth rate [mg/l = g/m³]
ρᵣ  = 6.0 * 1000 # mRNA density [g/l -> g/m³]
σ   = 0.345 * 86400# Biofilm detachment coefficient [1/m⋅s -> 1/m⋅d]
d   = 10 # Turnover rate of glucose below penetration depth Rp

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Acid Stress Response",
    tFinal=100*t₀,        # Simulation time [days]
    tol=1e-4,        # Tolerance
    outPeriod=0.1,   # Time between outputs [days]
    plotSize=(900,600), # Plot size [pixels]

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames=["Acid Stress Gene"],     # Particulate names
    Xto=[X],                    # Tank particulate concentration initial condition(s)
    Pbo=[X/ρᵣ], #???          # Biofilm particulates volume fraction initial condition(s) 
    rho=[ρᵣ],                 # Particulate densities
    Kdet=σ,                   # Particulates detachment coefficient
    srcX=[(S,X,Lf,t,z,p) -> -m*X], # Source of particulates
    # Growthrates for each particulate 
    #                                            growth when p<3000       else 0.0
    mu=[(S,X,Lf,t,z,p) -> map(i -> S[2,i] < μₚ ? μᵢ*S[1,i]*(1-S[2,i]/μₚ) : 0.0 ,1:length(S[1,:])) ],
    #mu=[(S,X,Lf,t,z,p) -> μᵢ*S[1,:].*(1.0 .- S[2,:]./μₚ)],

    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    SNames=["Glucose", "Lactate"], # Substrate names
    Sin=[(t) -> s⁰, (t) -> 0.0],   # Substrate inflow (can be function of time)
    Sto=[s⁰,0.0],            # Tank substrate concentration initial condition(s)
    Sbo=[0.0,0.0],           # Biofilm substrates concentration initial condition(s)
    Yxs=[Yxs -Yxs/Yps],      # Biomass yield coefficient on substrate
    Daq=[Ds, Dp],            # Substrate diffusion through boundary layer
    De =[Ds, Dp],            # Substrate diffusion through biofilm     
     # Source of substrates
    # removal of glucose when z < Lf-Rp, else 0.0
    srcS=[(S,X,Lf,t,z,p) -> map(i -> z[i] < Lf - Rp ? -d*S[1,i] : 0.0 ,1:length(S[1,:])), 
          (S,X,Lf,t,z,p) -> 0.0], 
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.1,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=1,          # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=50,          # Number of grid points in biofilm
    Lfo=L₀,     # Biofilm initial thickness [m]
    LL=1.0e-4,      # Boundary layer thickness [m]
)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,Xt,St,Pb,Sb,Lf,p) # Plot final results

savefig("Case2.pdf")