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

# Create empty dictionary to hold parameters 
d = createDict()

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "Biomass-Glucose-Lactate")
addParam!(d, "tFinal",   10)    # Simulation time [days]
addParam!(d, "tol",      1e-4)  # Tolerance
addParam!(d, "outPeriod",1)     # Time between outputs [days]
addParam!(d, "plotSize",(900,600)) # Plot size [pixels]

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Biomass"])    # Particulate names
addParam!(d, "Xto",   [Xt])           # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [1/6])          # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [ρᵣ])           # Particulate densities
addParam!(d, "Kdet",  σ)              # Particulates detachment coefficient
addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> 0.0])     # Source of particulates
# Growthrate: call μ(s,p) for S = [s,p]
addParam!(d, "mu", [(S,X,Lf,t,z,p) -> μ(S[1],S[2])])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Glucose", "Lactate"])     # Substrate names
addParam!(d, "Sin",   [(t) -> s⁰, (t) -> 0.0])    # Substrate inflow (can be function of time)
addParam!(d, "Sto",   [ s⁰, 0.0])      # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [0.0, 0.0])      # Biofilm substrates concentration initial condition(s)
addParam!(d, "Yxs",   [Yxs -Yxs/Yps])  # Biomass yield coefficient on substrate
addParam!(d, "Dt",    [ Ds, Dp])       # Aquious substrate diffusion through tank fluid
addParam!(d, "Db",    [ Ds, Dp])       # Effective substrate diffusion through biofilm

addParam!(d, "srcS",  [(S,X,Lf,t,z,p) -> 0.0,     # Source of substrates
                       (S,X,Lf,t,z,p) -> 0.0 ]) 

# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d, "V", 0.1)        # Volume of tank [m³]
addParam!(d, "A",   1)        # Surface area of biofilm [m²]
addParam!(d, "Q",   2)        # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d, "Nz",  50)      # Number of grid points in biofilm
addParam!(d, "Lfo", L₀)      # Biofilm initial thickness [m]
addParam!(d, "LL",  1.0E-4)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

# Run with inhibition 
####################
p_max = 400
t_in,zm_in,Xt_in,St_in,Pb_in,Sb_in,Lf_in,sol_in = BiofilmSolver(p)
plt = biofilm_plot(sol_in,p,"Inhibition",size=(900,600), line=(:blue,[:solid :dash]))

# Run without inhibition 
####################
p_max  = Inf  
t_no,zm_no,Xt_no,St_no,Pb_no,Sb_no,Lf_no,sol_no = BiofilmSolver(p)
# Add to plot
plt = biofilm_plot!(plt,sol_no,p,"No Inhibition",size=(900,600), line=(:red,[:solid :dash]))
display(plt)
savefig("Case2.pdf")

