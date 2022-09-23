using Biofilm 

# Constants used for growthrates of particulate(s)
KmB1 = 0.200; KmB3 = 11;  KmC2 = 20;  KI = 1;
mumaxA = 0.4;  mumaxB = 6.72;  mumaxC = 1.46;

# Define a structure to hold all the parameters
p = param(
    # Growthrates for each particulate (constants defined above!)
    mu=[(S, X, Lf, t, z, p) -> mumaxB*(S[1,:]./(KmB1.+S[1,:])).*(S[3,:]./(KmB3.+S[3,:])), # SOB
        (S, X, Lf, t, z, p) -> mumaxC*(S[2,:]./(KmC2.+S[2,:])).*(1.0./(1.0.+S[1,:]/KI))], # SRB

    # Source of particulates (constants defined above)
    src=[(S, X, t, p) -> 0.0,
         (S, X, t, p) -> 0.0,
         (S, X, t, p) -> 0.0],

    # Substrate inflow (can be function of time)
    Sin=[(t) -> 8.6
         (t) -> 48.0
         (t) -> 0.0],

    # Time
    tFinal=500,   # Simulation time [days]
    outPeriod=1,  # Time between outputs [days]
    plotPeriod=5,    # Time between plots [days] (make multiple of outPeriod!)

    # Simulation
    Title="SOB-SRB Test Case",
    SNames=["Oxygen","Sulfate","Hydrogen Sulfide"],
    XNames=["SOB - Sulfide-Oxidizer","SRB - Sulfate-Reducer"],
    #                           SOB                 SRB
    makePlots=true,

    # Tank Geometry
    V=0.01,        # Volume of tank [m³]
    A=1,          # Surface area of biofilm [m²]
    Q=10,          # Flowrate through tank [m³/s]
    Xo=[1.0e-6,1.0e-6],# Tank particulate initial condition(s)
    So=[8.6,48.0,0.0],    # Tank substrate initial condition(s)
    LL=2.0e-4,    # Boundary layer thickness [m]

    # Biofilm
    Nz=20,            # Number of grid points to represent biofilm
    Pbo=[0.05,0.15],     # Biofilm particulates initial condition(s)
    Sbo=[8.6,48.0,1e-5],     # Biofilm substrates initial condition(s)
    Lfo=5.0E-6,    # Biofilm initial thickness [m]

    # Substance Constants
    # Biomass yield coeffficient on substrate
    #     oxygen  sulfate  Hy. sulfide
    Yxs=[  0.58    0.0      0.09       # SOB uses oxygen and sulfide
           0.00    0.584   -1.645],    # SRB uses sulfate and produces sulfide
    Daq=[1.51e-4,8e-5,1.21e-4],    # Substrate diffusion through boundary layer
    De =[6.8e-5,4e-5,6.04e-5],     # Substrate diffusion through biofilm     
    rho=[2.5e5,2.5e5],     # Particulate densities
    Kdet=30.0,     # Particulates detachment coefficient

    # Tolerance
    tol=1e-4,
)

t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
makePlots(t,zm,X,S,Pb,Sb,Lf,p) # Plot final results

# Post processing
# Call movieBiofilm() to make a movie of the biofilm particulates and substrates
using Plots
function movieBiofilm()
    
    times=1:50 ### Adjust the times as needed ####

    # Make animation
    anim = @animate for t in times 
        analyzeBiofilm(sol,p,t,makePlot=true)
    end

    # View annimation 
    gif(anim) # Make a gif
    #gif(anim,"anim.mp4")  # Make a .mp4

end
#movieBiofilm()