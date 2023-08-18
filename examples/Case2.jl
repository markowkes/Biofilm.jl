using Biofilm 
using Plots
using LaTeXStrings
using Printf
using Measures

# This case is adapted from the Acid Stress Response simulation in https://doi.org/10.1371/journal.pone.0083626

# Notation
# s = Glucose
# p = Lactate

# Constants from paper (with unit conversions)
Ds  = 2.24e-10 * 86400 # Glucose diffusion coefficient [m²/s -> m²/d]
Dp  = 2.24e-10 * 86400 # Lactate diffusion coefficient [m²/s -> m²/d]
L₀  = 80.0 * 1e-6 # Initial biofilm thickness [μm -> mm]
Xt  = 5.0  # Biomass initial cell density [g/m³]
Yxs = 0.5 # Yield coefficient of biomass on glucose g_X/g_s
Yps = 0.9 # Yield coefficient of lactate on glucose g_p/g_s 
s⁰  = 800  # Bulk concentration of glucose [mg/l = g/m³]
μ_max  = 0.03  # Specific growth rate coefficient [1/d⋅m³/g]
p_max = 400.0 # Inhibition cutoff used in growth rate [mg/l = g/m³]
ρᵣ  = 3e5 # Biomass density [g/l -> g/m³]
# Adjust σ to get Lf[end] ≈ 400 μm
σ   = 1500 # Biofilm detachment coefficient [1/m⋅s -> 1/m⋅d]

# Define mu function 
function μ(s,p)
    if p < p_max 
        μ = μ_max*s*(1-p/p_max)
    else
        μ = 0.0
    end
    return μ
end

# Tuple to hold all input parameters
p = (
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title =     "Biomass-Glucose-Lactate",
    tFinal =    10,    # Simulation time [days]
    tol =       1e-4,  # Tolerance
    outPeriod = 1,     # Time between outputs [days]
    plotSize = (900,600), # Plot size [pixels]

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    XNames = ["Biomass"],    # Particulate names
    Xto =    [Xt],           # Tank particulate concentration initial condition(s)
    Pbo =    [1/6],          # Biofilm particulates volume fraction initial condition(s) 
    rho =    [ρᵣ],           # Particulate densities
    Kdet =   σ,              # Particulates detachment coefficient
    srcX =   [(S,X,Lf,t,z,p) -> 0.0],     # Source of particulates
    # Growthrate: call μ(s,p) for S = [s,p]
    mu = [(S,X,Lf,t,z,p) -> μ(S[1],S[2])],

    # ----------------- #
    # Solute Parameters #
    # ----------------- #
    SNames = ["Glucose", "Lactate"],   # Solute names
    Sin =  [(t) -> s⁰, (t) -> 0.0],    # Solute inflow (can be function of time)
    Sto =  [ s⁰, 0.0],      # Tank substrate concentration initial condition(s)
    Sbo =  [0.0, 0.0],      # Biofilm substrates concentration initial condition(s)
    Yxs =  [Yxs -Yxs/Yps],  # Biomass yield coefficient on substrate
    Dt =   [ Ds, Dp],       # Aquious substrate diffusion through tank fluid
    Db =   [ Ds, Dp],       # Effective substrate diffusion through biofilm
    srcS = [(S,X,Lf,t,z,p) -> 0.0,
            (S,X,Lf,t,z,p) -> 0.0 ],

    # --------------- #
    # Tank Parameters #
    # --------------- #
    V = 0.1,        # Volume of tank [m³]
    A =   1,        # Surface area of biofilm [m²]
    Q =   2,        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz =  50,      # Number of grid points in biofilm
    Lfo = L₀,      # Biofilm initial thickness [m]
    LL =  1.0E-4,  # Boundary layer thickness [m]
)

# Run with inhibition 
####################
p_max = 400
t_in,zm_in,Xt_in,St_in,Pb_in,Sb_in,Lf_in,sol_in = BiofilmSolver(p)
plt = biofilm_plot(sol_in,p,"Inhibition",size=(900,600), line=(:blue,[:solid :dashdot]))

# Run without inhibition 
####################
p_max  = Inf  
t_no,zm_no,Xt_no,St_no,Pb_no,Sb_no,Lf_no,sol_no = BiofilmSolver(p)
# Add to plot
plt = biofilm_plot!(plt,sol_no,p,"No Inhibition",size=(900,600), line=(:red,[:solid :dashdot]))
display(plt)
savefig("Case2.pdf")

