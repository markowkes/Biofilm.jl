
using Biofilm
using UnPack

# Constants used for growthrates of particulate(s)
mumax = 20;
KM = 3;

# Define a structure to hold all the parameters
p = param(
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> (mumax * S[1,:]) ./ (KM .+ S[1,:])],

    # Source of particulates
    src=[(S, X, p) -> 0.0],

    # Substrate inflow (can be function of time)
    Sin=[(t) -> 25, (t) -> 25],

    # Time
    tFinal=30,   # Simulation time [days]
    outPeriod=1e-1,  # Time between outputs [days]

    # Simulation
    Title="Multiple Independent S Case",
    SNames=["Substrate 1","Substrate 2"],
    XNames=["Bug"],
    makePlots=true,

    # Tank Geometry
    V=0.1,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=1,          # Flowrate through tank [m³/s]
    Xo=[10.0],     # Tank particulate initial condition(s)
    So=[25.0,25.0],    # Tank substrate initial condition(s)
    LL=1.0e-4,    # Boundary layer thickness [m]

    # Biofilm
    Nz=50,            # Number of grid points to represent biofilm
    Pbo=[0.2],     # Biofilm particulates initial condition(s)
    Sbo=[0.0,0.0],     # Biofilm substrates initial condition(s)
    Lfo=5.0E-6,    # Biofilm initial thickness [m]

    # Substance Constants
    Yxs=[0.5 Inf]::Array{Float64,2},     # Biomass yield coeffficient on substrate
    Daq=[4.0e-5, 6.0e-5],    # Substrate diffusion through boundary layer
    De =[1.0e-5, 1.5e-5],    # Substrate diffusion through biofilm     
    rho=[1.0e5],     # Particulate densities
    Kdet=1900.0,     # Particulates detachment coefficient

    # Tolerance
    tol=1e-2,
)

@unpack Nx,Ns,Nz = p

# Create test solution
X=p.Xo
S=p.So
Pb=p.Pbo.*ones(Nx,Nz)
Xb=Pb./p.rho
Sb=p.Sbo.*ones(Ns,Nz)
Lf=p.Lfo
t=0.5

# Create test grid
z=range(0.0,Lf,Nz+1)
zm=0.5*(z[1:Nz]+z[2:Nz+1])
dz=z[2]-z[1]
g=Biofilm.biofilmGrid(z,zm,dz)

# Compute Flux S
fluxS = Biofilm.computeFluxS(S,Sb,p,g)              # Flux of substrate in biofilm
μb    = Biofilm.computeMu_biofilm(Sb,Xb,Lf,t,p,g)   # Growthrates in biofilm
μt    = Biofilm.computeMu_tank(S,X,Lf,t,p,g)        # Growthrates in tank
V     = Biofilm.computeVel(μb,Pb,p,g)               # Velocity of particulates
Vdet  = p.Kdet*Lf^2                                 # Detachment velocity
fluxP = Biofilm.computeFluxP(Pb,V,Vdet,p)           # Flux of particulates in biofilm

dX    = Biofilm.dXdt(t,X,S,Xb,Lf,Vdet,μt,p,g)
dS    = Biofilm.dSdt(t,X,S,Lf,Pb,fluxS,μt,p)        # Tank substrates
dPb   = Biofilm.dPbdt(μb,Sb,Pb,fluxP,p,g)           # Biofilm particulates 
dSb   = Biofilm.dSbdt(μb,Xb,fluxS,p,g)              # Biofilm substrates
dLf   = Biofilm.dLfdt(V,Vdet)                       # Biofilm thickness

