# Create figures in parts for presentation overlays

begin # Load packages 
    using Biofilm 
    using OrdinaryDiffEq
    using Measures
    using Plots
    using Printf
    using Statistics
    using UnPack
    using JLD2
    using Accessors
end

begin # Setup biofilm

    # Create empty dictionary to hold parameters 
    d = createDict()

    smoothHeaviside(t,t0)=0.5*tanh.(10*(t.-t0).-0.5).+0.5

    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    addParam!(d, "Title",    "Hydrogen Peroxide Dosing")
    addParam!(d, "tFinal",  100)    # Simulation time [days]
    addParam!(d, "tol",      1e-8)  # Tolerance
    addParam!(d, "outPeriod",10.0)   # Time between outputs [days]
    addParam!(d, "plotPeriod",10)   # Time between plots [days] 
    #addParam!(d, "discontinuityPeriod",2.5) # Let solver know when discontinuities
    addParam!(d, "makePlots",false)

    # ---------------------- #
    # Particulate Parameters #
    # ---------------------- #
    addParam!(d, "XNames",["Live","Dead"])    # Particulate names
    addParam!(d, "Xto",   [1.0, 0.0])         # Tank particulate concentration initial condition(s)
    addParam!(d, "Pbo",   [0.08, 0.0])        # Biofilm particulates volume fraction initial condition(s) 
    addParam!(d, "rho",   [2.5e5, 2.5e5])     # Particulate densities
    addParam!(d, "Kdet",  1e4)                # Particulates detachment coefficient
    k_dis = 0.5; #m³/g/d   # Source term constant
    addParam!(d, "srcX",  [(S,X,Lf,t,z,p) -> -k_dis*X[1].*S[2], # Source of particulates
                        (S,X,Lf,t,z,p) -> +k_dis*X[1].*S[2] ])
    # Growthrates for each particulate
    mumax = 9.6; # 1/days
    KM = 5; # g/m³
    addParam!(d, "mu", [(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM .+ S[1]) 
                        (S,X,Lf,t,z,p) -> 0.0 ] )

    # -------------------- #
    # Substrate Parameters #
    # -------------------- #
    addParam!(d, "SNames",["Glucose","Hydrogen Peroxide"])     # Substrate names
    GlucoseIn = 100; #g/m³
    dose1 = 500.0; dose2 = 0.0;
    addParam!(d, "Sin",   [
            (t) -> GlucoseIn,        # Substrate inflow (can be function of time)
            (t) -> dose1*smoothHeaviside(t,2.0)+dose2*smoothHeaviside(t,6.0)])
    addParam!(d, "Sto",   [100.0, 0.0])          # Tank substrate concentration initial condition(s)
    addParam!(d, "Sbo",   [  0.0, 0.0])          # Biofilm substrates concentration initial condition(s)
    # Biomass yield coefficient on substrate
    addParam!(d, "Yxs",   [#Glucose   H. Per.
                            0.26       0.0       # Live use glucose
                            0.00       0.0   ])  # Dead doesn't use/produce anything
    addParam!(d, "Daq",   [5.2e-5, 1.09e-4])     # Substrate diffusion through boundary layer
    addParam!(d, "De",    [1.3e-5, 6.52e-5])     # Substrate diffusion through biofilm     
    k_bl  = 10.0; #m³/g/d
    k_bd  = 10.0; #m³/g/d
    addParam!(d, "srcS",  [                      # Source of substrates
        (S,X,Lf,t,z,p) -> 0.0,  
        (S,X,Lf,t,z,p) -> -k_bl*X[1].*S[2].-k_bd*X[2].*S[2] ])

    # --------------- #
    # Tank Parameters #
    # --------------- #
    addParam!(d, "V", 0.1)        # Volume of tank [m³]
    addParam!(d, "A",   1)        # Surface area of biofilm [m²]
    addParam!(d, "Q",   1)        # Flowrate through tank [m³/d]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    addParam!(d, "Nz",  50)        # Number of grid points in biofilm
    addParam!(d, "Lfo", 50.0e-6)   # Biofilm initial thickness [m]
    addParam!(d, "LL",  1.0e-5)    # Boundary layer thickness [m]

    # Package and check parameters 
    p = packageCheckParam(d)

end

#########
# Fig 1: Biofilm thickness dynamics. A, no biocide treatment. B, biocide dosed continuously beginning at t = 2 days. C, biocide dosed for 4 days beginning at t = 2 days, then removed to allow biofilm recovery.
#########
# Run 3 cases
begin # runs 
    # Run cases if needed
    if !(@isdefined(sol1) && @isdefined(sol2) && @isdefined(sol3))
        k_bl  = 10.0; #m³/g/d
        k_bd  = 10.0; #m³/g/d
        GlucoseIn = 100.0; # g/m³
        # Run out to 2000 days for Fig. 4a 
        # (parameter studies will only run to 100 days)
        p2000 = @set p.tFinal = 2000
        dose1 =   0.0; dose2 =    0.0; t1,zm1,X1,S1,Pb1,Sb1,Lf1,sol1 = BiofilmSolver(p2000)
        dose1 = 500.0; dose2 =    0.0; t2,zm2,X2,S2,Pb2,Sb2,Lf2,sol2 = BiofilmSolver(p2000)
        dose1 = 500.0; dose2 = -500.0; t3,zm3,X3,S3,Pb3,Sb3,Lf3,sol3 = BiofilmSolver(p2000)
    end
end
begin # plots
    # Part 1
    fig = plot()
    fig = plot!(
            xlabel = "Time (d)",
            ylabel = "Biofilm Thickness (μm)",
            ylims = (0.0,200),
            xlims = (0,10),
            legend=:bottomright,
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            foreground_color_legend = nothing,
            size=(650,400),
        )
    fig = plot!(t1,Lf1'.*1e6,linecolor = :black,line=(:solid,2), label="No Dosing     ")
    # Blank lines
    fig = plot!(t1,-Lf1'.*1e6,linecolor = :white, label=" ")
    fig = plot!(t1,-Lf1'.*1e6,linecolor = :white, label=" ")
    savefig("Fig1_part1.pdf")
    
    # Part 2
    fig = plot()
    fig = plot!(
            xlabel = "Time (d)",
            ylabel = "Biofilm Thickness (μm)",
            ylims = (0.0,200),
            xlims = (0,10),
            legend=:bottomright,
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            foreground_color_legend = nothing,
            size=(650,400),
        )
    fig = plot!(t1,Lf1'.*1e6,linecolor = :black,line=(:solid,2), label="No Dosing     ")
    fig = plot!(t2,Lf2'.*1e6,linecolor = :blue ,line=(:dash, 2), label="Dosing on")
    fig = quiver!([2], [100], quiver=([0], [30]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([2],[90],"Dosing On")
    # Blank lines
    fig = plot!(t1,-Lf1'.*1e6,linecolor = :white, label=" ")
    savefig("Fig1_part2.pdf")

    # Part 3
    fig = plot()
    fig = plot!(
            xlabel = "Time (d)",
            ylabel = "Biofilm Thickness (μm)",
            ylims = (0.0,200),
            xlims = (0,10),
            legend=:bottomright,
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            foreground_color_legend = nothing,
            size=(650,400),
        )
    fig = plot!(t1,Lf1'.*1e6,linecolor = :black,line=(:solid,2), label="No Dosing     ")
    fig = plot!(t2,Lf2'.*1e6,linecolor = :blue ,line=(:dash, 2), label="Dosing on")
    fig = quiver!([2], [100], quiver=([0], [30]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([2],[90],"Dosing On")
    fig = plot!(t3,Lf3'.*1e6,linecolor = :red  ,line=(:dot,  2), label="Dosing on/off")
    fig = quiver!([6], [130], quiver=([0], [-30]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([6],[90],"Dosing Off")
    savefig("Fig1_part3.pdf")
    display(fig)
end


#########
# Fig 5: Biofilm protection depends on hydrogen peroxide neutralizing activity of dead cells. A, biofilm thickness versus time with neutralization of hydrogen peroxide turned off for either live or dead cells. B, percent live cells versus time with neutralization of hydrogen peroxide turned off for either live or dead cells.
#This is relative to base case
#########
begin # run
    if !(@isdefined(sol4) && @isdefined(sol5))
        GlucoseIn = 100.0; # g/m³

        # Turn off live neutralization (only dead)
        dose1 = 500.0; dose2 = 0.0; k_bl=0.0; k_bd=10.0;
        t4,zm4,X4,S4,Pb4,Sb4,Lf4,sol4 = BiofilmSolver(p)

        # Turn off dead neutralization (only live)
        dose1 = 500.0; dose2 = 0.0; k_bl=10.0; k_bd=0.0;
        t5,zm5,X5,S5,Pb5,Sb5,Lf5,sol5 = BiofilmSolver(p)
    end
end
begin # plot Biofilm Thickness vs Time
    #Part 1
    fig = plot()
    fig = plot!(
        xlabel = "Time (d)",
        ylabel = "Biofilm Thickness (μm)",
        ylims = (0.0,200),
        xlims = (0,10),
        legend=:right,
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (800,500),
        margin = 10mm,
    )
    fig = plot!(t2,Lf2'.*1e6,line=(:solid,2, :black), label="Live & Dead Neutralization")
    fig = quiver!([2], [20], quiver=([0], [40]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([1.9],[40],text("Dosing On",:right,12))
    # Add blank lines for legend
    fig = plot!(t2,-Lf2'.*1e6,line=(:white),label=" ")
    fig = plot!(t2,-Lf2'.*1e6,line=(:white),label="  ")
    savefig("Fig5a_part1.pdf")


    #Part 2
    fig = plot()
    fig = plot!(
        xlabel = "Time (d)",
        ylabel = "Biofilm Thickness (μm)",
        ylims = (0.0,200),
        xlims = (0,10),
        legend=:right,
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (800,500),
        margin = 10mm,
    )
    fig = plot!(t2,Lf2'.*1e6,line=(:solid,2, :black), label="Live & Dead Neutralization")
    fig = plot!(t4,Lf4'.*1e6,line=(:dash, 2, :purple),label="Only Dead Neutralization")
    fig = quiver!([2], [20], quiver=([0], [40]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([1.9],[40],text("Dosing On",:right,12))
    # Add blank lines for legend
    fig = plot!(t2,-Lf2'.*1e6,line=(:white),label="   ")
    savefig("Fig5a_part2.pdf")

    #Part 3
    fig = plot()
    fig = plot!(
        xlabel = "Time (d)",
        ylabel = "Biofilm Thickness (μm)",
        ylims = (0.0,200),
        xlims = (0,10),
        legend=:right,
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (800,500),
        margin = 10mm,
    )
    fig = plot!(t2,Lf2'.*1e6,line=(:solid,2, :black), label="Live & Dead Neutralization")
    fig = plot!(t4,Lf4'.*1e6,line=(:dash, 2, :purple),label="Only Dead Neutralization")
    fig = plot!(t5,Lf5'.*1e6,line=(:dot,  2, :green), label="Only Live Neutralization")
    fig = quiver!([2], [20], quiver=([0], [40]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([1.9],[40],text("Dosing On",:right,12))
    savefig("Fig5a_part3.pdf")

end