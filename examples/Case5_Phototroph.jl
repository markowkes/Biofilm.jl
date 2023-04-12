using Biofilm 

# Create empty dictionary to hold parameters 
d = createDict()

# Define light as a function of time and depth within biofilm
diss=2000;  # Dissipation rate into biofilm [1/m]
smoothHeaviside(t,t0)=0.5*tanh.(100*(t.-t0).-0.5).+0.5
# Light :              turns off at t=0.25             turns on at t=0.75
intensity(t) = smoothHeaviside(mod(t,1),0.25)-smoothHeaviside(mod(t,1),0.75)
# Dissipation of light into biofilm (1 at top with a rate of decrease of diss)
dissipation(z,Lf) = max.(0.0,1.0.-(Lf.-z)*diss)
light(t,z,Lf) = intensity(t)*dissipation(z,Lf)

# Check light intensity with time 
using Plots
Lf = 600e-6; z = Lf;
t = range(0.0,3.0,1000)
plot(t,map(t -> light(t,z,Lf),t),
    xlabel="Time [days]",
    ylabel="Light Intensity",
    legend=false,
    size=(400,200),
    )
savefig("LightTime.pdf")

# Check light intensity with depth 
using Plots
Lf = 600e-6; t = 0.5; 
z = range(0.0, Lf, 1000)
plot(1e6.*z,map(z -> light(t,z,Lf),z),
    xlabel="Height in Biofilm [μm]",
    ylabel="Light Intensity",
    legend=false,
    size=(400,200),
    )
savefig("LightDepth.pdf")

# --------------------- #
# Simulation Parameters #
# --------------------- #
addParam!(d, "Title",    "Phototroph Case")
addParam!(d, "tFinal",   50)   # Simulation time [days]
addParam!(d, "tol",      1e-4)  # Tolerance
addParam!(d, "outPeriod",5.0)   # Time between outputs [days]
# Let solver know when discontinuities (changes in light) occur
addParam!(d, "discontinuityPeriod",0.25)  

# ---------------------- #
# Particulate Parameters #
# ---------------------- #
addParam!(d, "XNames",["Phototroph"])    # Particulate names
addParam!(d, "Xto",   [1.0])             # Tank particulate concentration initial condition(s)
addParam!(d, "Pbo",   [0.2])             # Biofilm particulates volume fraction initial condition(s) 
addParam!(d, "rho",   [2.5e5])           # Particulate densities
addParam!(d, "Kdet",  100.0)             # Particulates detachment coefficient
b = 0.1 # Source term constant
addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> 0.0]) # Source of particulates
# Growthrates for each particulate
# Constants used for growthrates of particulate(s)
mumax = 0.4;
addParam!(d, "mu", [(S,X,Lf,t,z,p) -> mumax*light(t,z,Lf)])

# -------------------- #
# Substrate Parameters #
# -------------------- #
addParam!(d, "SNames",["Oxygen"])     # Substrate names
addParam!(d, "Sin",   [(t) -> 8.6])   # Substrate inflow (can be function of time)
addParam!(d, "Sto",   [8.6])          # Tank substrate concentration initial condition(s)
addParam!(d, "Sbo",   [8.6])          # Biofilm substrates concentration initial condition(s)
addParam!(d, "Yxs",   [-0.52])        # Biomass yield coefficient on substrate
addParam!(d, "Dt",    [1.51e-4])      # Aquious substrate diffusion through tank fluid
addParam!(d, "Db",    [6.8E-5])       # Effective substrate diffusion through biofilm
addParam!(d, "srcS",  [(S,X,Lf,t,z,p) -> 0.0])     # Source of substrates

# --------------- #
# Tank Parameters #
# --------------- #
addParam!(d, "V", 0.01)       # Volume of tank [m³]
addParam!(d, "A",   1)        # Surface area of biofilm [m²]
addParam!(d, "Q",  10)        # Flowrate through tank [m³/d]

# ------------------ #
# Biofilm Parameters #
# ------------------ #
addParam!(d, "Nz",  50)       # Number of grid points in biofilm
addParam!(d, "Lfo", 5.0e-6)   # Biofilm initial thickness [m]
addParam!(d, "LL",  2.0e-4)  # Boundary layer thickness [m]

# Package and check parameters 
p = packageCheckParam(d)

t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver
biofilm_plot(sol,p,size=(900,600))
savefig("Case5.pdf")

# Postprocessing Results 
# ----------------------

# Biofilm quantities 
tout = 49.5 # Time when light is on
biofilm_analyze(sol,p,tout,makePlot=true,plotSize=(900,325))
# Could also make plot by directly calling recipe 
# biofilm_plot_film(sol([0,tout]),p,size=(900,325))
savefig("Case5_lighton.pdf")

# Times to analyze solution 
tout = 48.0:0.01:50.0; 
# Get solution at certain times
Xtout,Stout,Lfout = biofilm_analyze(sol,p,tout)

# Function to plot when the light is on
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
function plot_light()
    plot(rectangle(0.5,1000.0,48.25,0.0), 
        opacity=.2, 
        fillcolor=:yellow,
        linecolor=:yellow,
        )
    plot!(rectangle(0.5,1000.0,49.25,0.0), 
        opacity=.2, 
        fillcolor=:yellow,
        linecolor=:yellow,
        )
end

# Plot tank particulate concentration versus time 
plot_light()
plot!(tout,Xtout',
    linecolor=:blue,
    xlabel=("Time [days]"),
    ylabel=("Phototroph Conc. [g/m³]"),
    legend=false,
    size=(300,300),
    ylims=(minimum(Xtout)-0.001,maximum(Xtout)+0.001),
    )
savefig("Case5_Xt.pdf")

# Plot tank substrate concentration versus time 
plot_light()
plot!(tout,Stout',
    linecolor=:blue,
    xlabel=("Time [days]"),
    ylabel=("Oxygen Conc. [g/m³]"),
    legend=false,
    size=(300,300),
    ylims=(minimum(Stout)-0.1,maximum(Stout)+0.1),
    )
savefig("Case5_St.pdf")

# Plot biofilm thickness versus time 
plot_light()
plot!(tout,1e6.*Lfout,
    linecolor=:blue,
    xlabel=("Time [days]"),
    ylabel=("Biofilm Thickness [μm]"),
    legend=false,
    size=(300,300),
    ylims=(minimum(1e6*Lfout)-2,maximum(1e6*Lfout)+2),
    )
savefig("Case5_Lf.pdf")


