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
addParam!(d, "Sto",   [ s⁰, 0.0])       # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [0.0, 0.0])       # Biofilm substrates concentration initial condition(s)
addParam!(d, "Yxs",   [Yxs -Yxs/Yps])  # Biomass yield coefficient on substrate
addParam!(d, "Daq",   [ Ds, Dp])        # Substrate diffusion through boundary layer
addParam!(d, "De",    [ Ds, Dp])        # Substrate diffusion through biofilm     
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

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver

# Run with inhibition 
####################
t_in,zm_in,Xt_in,St_in,Pb_in,Sb_in,Lf_in,sol_in = BiofilmSolver(p)
savefig("Case2.pdf")

# Compute growthrates for plots
Xb_in=similar(Pb_in)
for j=1:p.Nx
    Xb_in[j,:] = p.rho[j]*Pb_in[j,:]  # Compute particulate concentrations
end
g_in = Biofilm.computeGrid(Lf_in[end],p)
μb_in = Biofilm.computeMu_biofilm(Sb_in,Xb_in,Lf_in[end],t_in[end],p,g_in)

# Run without inhibition 
####################
p_max  = Inf  
t_no,zm_no,Xt_no,St_no,Pb_no,Sb_no,Lf_no,sol_no = BiofilmSolver(p)
# Compute growthrates for plots
Xb_no=similar(Pb_no)
for j=1:p.Nx
    Xb_no[j,:] = p.rho[j]*Pb_no[j,:]  # Compute particulate concentrations
end
g_no = Biofilm.computeGrid(Lf_no[end],p)
μb_no = Biofilm.computeMu_biofilm(Sb_no,Xb_no,Lf_no[end],t_no[end],p,g_no)

# Make comparison plots 
# (similar to makePlots in outputs.jl)
#############################################

# Adjust names to work with legends
p.Nx==1 ? Xs=p.XNames[1] : Xs=reshape(p.XNames,1,length(p.XNames))
p.Ns==1 ? Ss=p.SNames[1] : Ss=reshape(p.SNames,1,length(p.SNames))

# Tank particulate concentration
p1= plot(t_in,Xt_in',label=Xs.*" : Inhibition",   line=:blue)
p1=plot!(t_no,Xt_no',label=Xs.*" : No Inhibition",line=:red)
xaxis!(L"\textrm{Time~[days]}")
yaxis!(L"\textrm{Tank~Particulate~Conc.~} [g/m^3]")

# Tank substrate concentration
p2= plot(t_in,St_in',label=Ss.*" : Inhibition",   line=(:blue,[:solid :dash]))
p2=plot!(t_no,St_no',label=Ss.*" : No Inhibition",line=( :red,[:solid :dash]))
xaxis!(L"\textrm{Time~[days]}")
yaxis!(L"\textrm{Tank~Substrate~Conc.~} [g/m^3]")

# Biofilm thickness
p3= plot(t_in,Lf_in'*1e6,label="Thickness : Inhibition",   line=:blue)
p3=plot!(t_no,Lf_no'*1e6,label="Thickness : No Inhibition",line=:red)
xaxis!(L"\textrm{Time~[days]}")
yaxis!(L"\textrm{Biofilm~Thickness~} [μm]")

# Biofilm particulate volume fractioin 
p4= plot(1e6.*zm_in,Pb_in',label=Xs.*" : Inhibition",   line=:blue,ylim=Biofilm.pad_ylim(Pb_in'))
p4=plot!(1e6.*zm_no,Pb_no',label=Xs.*" : No Inhibition",line=:red ,ylim=Biofilm.pad_ylim(Pb_no'))
#p4=plot!(zm,sum(Pb,dims=1)',label="Sum")
xaxis!(L"\textrm{Height~in~Biofilm~} [\mu m]")
yaxis!(L"\textrm{Biofilm~Particulate~Vol.~Frac.~[-]}")

# Biofilm substrate concentration
p5= plot(1e6.*zm_in,Sb_in',label=Ss.*" : Inhibition",   line=(:blue,[:solid :dash]))
p5=plot!(1e6.*zm_no,Sb_no',label=Ss.*" : No Inhibition",line=( :red,[:solid :dash]))
xaxis!(L"\textrm{Height~in~Biofilm~} [\mu m]")
yaxis!(L"\textrm{Biofilm~Substrate~Conc.~} [g/m^3]")

# Particulate growthrates vs depth
p6= plot(1e6.*zm_in,μb_in',label=Xs.*" : Inhibition",   line=(:blue))
p6=plot!(1e6.*zm_no,μb_no',label=Xs.*" : No Inhibition",line=:red)
xaxis!(L"\textrm{Height~in~Biofilm~} [\mu m]")
yaxis!(L"\textrm{Particulate~Growthrates~[-]}")

# Put plots together
myplt=plot(p1,p2,p3,p4,p5,p6,
    layout=(2,3),
    size=p.plotSize,
    plot_title=@sprintf("%s : t = %.2f",p.Title,t_in[end]),
    #plot_titlevspan=0.5,
    left_margin=10mm, 
    bottom_margin=10mm,
    foreground_color_legend = nothing,
    legend = :outertop,
)
display(myplt)

savefig("Case2_noinhibition.pdf")

