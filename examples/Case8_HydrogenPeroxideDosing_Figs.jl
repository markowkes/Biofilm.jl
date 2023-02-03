using Biofilm 
using DifferentialEquations
using Measures
using Plots
using Printf
using Statistics
using UnPack
using JLD2
using DifferentialEquations

# Constants used for growthrates of particulate(s)
Km = 5; # g/m³
mumax = 9.6; # 1/days
k_dis = 0.5; #m³/g/d
GlucoseIn = 100; #g/m³

smoothHeaviside(t,t0)=0.5*tanh.(1000*(t.-t0).-0.5).+0.5

# Define a structure to hold all the parameters
p = param(
    # --------------------- #
    # Simulation Parameters #
    # --------------------- #
    Title="Hydrogen Peroxide Dosing",
    tFinal=100,     # Simulation time [days]
    tol=1e-8,       # Tolerance
    outPeriod=2,    # Time between outputs [days]
    plotPeriod=10,    # Time between plots [days] (make multiple of outPeriod!)
    #optionalPlot="source", # 6th plot: "growthrate" (default) or "source"
    #discontinuityPeriod=2.5,  # Let solver know when discontinuities (changes in light) occur
    makePlots=false,

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
    Sin=[(t) -> GlucoseIn,    # Substrate inflow (can be function of time)
         (t) -> dose1*smoothHeaviside(t,2.0)+dose2*smoothHeaviside(t,6.0)],
    Sto=[100.0,0.0],  # Tank substrate concentration initial condition(s)
    Sbo=[0.0,0.0], # Biofilm substrates concentration initial condition(s)
    # Biomass yield coefficient on substrate
    #     Glucose   H. Per.
    Yxs=[ 0.26       0.0        # Live use glucose
          0.00       0.0   ],   # Dead doesn't use/produce anything
    Daq=[5.2e-5, 1.09e-4],    # Substrate diffusion through boundary layer
    De =[1.3e-5, 6.52e-5],     # Substrate diffusion through biofilm     
    srcS=[(S, X, t, p) -> 0.0,  
          (S, X, t, p) -> -k_bl*X[1,:].*S[2,:].-k_bd*X[2,:].*S[2,:] ],
           
    # --------------- #
    # Tank Parameters #
    # --------------- #
    V=0.1,        # Volume of tank [m³]
    A=1.0,        # Surface area of biofilm [m²]
    Q=1.0,        # Flowrate through tank [m³/s]

    # ------------------ #
    # Biofilm Parameters #
    # ------------------ #
    Nz=50,          # Number of grid points in biofilm
    Lfo=50.0e-6,     # Biofilm initial thickness [m]
    LL=1.0e-5,      # Boundary layer thickness [m]
)

"""
Compute the % Live averaged throughout biofilm 
"""
function meanLive(Pb)
    return mean(Pb[1,:])./(mean(Pb[1,:])+mean(Pb[2,:]))*100
end
function meanDead(Pb)
    return mean(Pb[2,:])./(mean(Pb[1,:])+mean(Pb[2,:]))*100
end


#########
# Fig 1: Biofilm thickness dynamics. A, no biocide treatment. B, biocide dosed continuously beginning at t = 2 days. C, biocide dosed for 4 days beginning at t = 2 days, then removed to allow biofilm recovery.
#########
# Run 3 cases
begin # runs 
    # Just run cases as sol is needed
    # if  isfile("cases1.jld2") && 
    #     isfile("cases2.jld2") && 
    #     isfile("cases3.jld2")
    #     println("Using saved files for cases1, cases2, and cases3")
    #     JLD2.@load "cases1.jld2" t1 zm1 X1 S1 Pb1 Sb1 Lf1 sol1
    #     JLD2.@load "cases2.jld2" t2 zm2 X2 S2 Pb2 Sb2 Lf2 sol2
    #     JLD2.@load "cases3.jld2" t3 zm3 X3 S3 Pb3 Sb3 Lf3 sol3
    # else
        k_bl  = 10.0; #m³/g/d
        k_bd  = 10.0; #m³/g/d
        GlucoseIn = 100.0; # g/m³
        dose1 =   0.0; dose2 =    0.0; t1,zm1,X1,S1,Pb1,Sb1,Lf1,sol1 = BiofilmSolver(p)
        dose1 = 500.0; dose2 =    0.0; t2,zm2,X2,S2,Pb2,Sb2,Lf2,sol2 = BiofilmSolver(p)
        dose1 = 500.0; dose2 = -500.0; t3,zm3,X3,S3,Pb3,Sb3,Lf3,sol3 = BiofilmSolver(p)
    #     JLD2.@save "cases1.jld2" t1 zm1 X1 S1 Pb1 Sb1 Lf1 sol1
    #     JLD2.@save "cases2.jld2" t2 zm2 X2 S2 Pb2 Sb2 Lf2 sol2
    #     JLD2.@save "cases3.jld2" t3 zm3 X3 S3 Pb3 Sb3 Lf3 sol3
    # end
end
begin # plots
    fig = plot()
    fig = plot!(t1,Lf1'.*1e6,linecolor = :black,line=(:solid,2), label="Dosing off")
    fig = plot!(t2,Lf2'.*1e6,linecolor = :blue ,line=(:dash, 2), label="Dosing on")
    fig = plot!(t3,Lf3'.*1e6,linecolor = :red  ,line=(:dot,  2), label="Dosing on/off")
    fig = quiver!([2], [100], quiver=([0], [30]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([2],[90],"Dosing On")
    fig = quiver!([6], [130], quiver=([0], [-30]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([6],[90],"Dosing Off")
    #fig = quiver!([2], [0], quiver=([0], [40]), linecolor = :black, line=(:solid, 1))
    # fig = annotate!([1.9],[20],text("Dosing On",:right,12))
    # fig = quiver!([6], [40], quiver=([0], [-40]), linecolor = :black, line=(:solid, 1))
    # fig = annotate!([5.9],[20],text("Dosing Off",:right,12))
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
        )
    display(fig)
    savefig("Fig1.svg")
    using Printf 
    @printf("Final biofilm thickness = %10.3f μm\n", Lf1[end]*1e6)
    @printf("Final biofilm thickness = %10.3f μm\n", Lf2[end]*1e6)
    @printf("Final biofilm thickness = %10.3f μm\n", Lf3[end]*1e6)
end

#########
# Fig 2: Steady state biofilm thickness, A, and percent live cells, B, when subjected to continuous biocide treatment.
#A, Plot Lf versus Cdose for concentrations ranging from about 1 g/m3 to about 20,000 g/m3.
#B, Plot % Live cells versus Cdose for same concentration range.
#########
begin # runs
    doses = 0.0:50.0:20000.0
    #doses = 50.0:5000.0:20000.0 # testing
    if  isfile("doses.jld2")
        println("Using saved files for doses")
        JLD2.@load "doses.jld2" ts zms Xs Ss Pbs Sbs Lfs
    else
        k_bl  = 10.0; #m³/g/d
        k_bd  = 10.0; #m³/g/d
        GlucoseIn = 100.0; # g/m³
        ts  = Vector{Float64}[]
        zms = Vector{Float64}[]
        Xs  = Matrix{Float64}[]
        Ss  = Matrix{Float64}[]
        Pbs = Matrix{Float64}[]
        Sbs = Matrix{Float64}[]
        Lfs = Matrix{Float64}[]
        for i in eachindex(doses)
            global dose1 = doses[i]
            global dose2 = 0.0
            println("Running solver with dose1=$dose1")
            t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p)
            push!( ts, t)
            push!(zms,collect(zm))
            push!( Xs, X)
            push!( Ss, S)
            push!(Pbs,Pb)
            push!(Sbs,Sb)
            push!(Lfs,Lf)
        end
        JLD2.@save "doses.jld2" ts zms Xs Ss Pbs Sbs Lfs
    end
end
begin # Plot Lf vs doses (add additional doses - make sure at S.S.)
    fig = plot()
    fig = plot!([0,maximum(doses)],[Lf1[end],Lf1[end]].*1e6, 
        linewidth = 1, 
        linecolor = :black,
        line=(:dot),
    )
    fig = annotate!(1.5e4,150,"No Dose Thickness")
    fig = plot!(doses,map(i -> Lfs[i][end]*1e6,1:length(doses)),
        linewidth = 2,
        linecolor = :black,
    )
    # Add label at crossover point 
    A = map(i -> Lfs[i][end]*1e6,1:length(doses))
    A = abs.(A .- A[1])
    A[1]=100.0
    value,index = findmin(A)
    xdose = doses[index]
    xLf   = Lfs[index][end]*1e6
    fig = annotate!([xdose],[xLf-50],@sprintf("%4.0f g/m³",xdose))
    fig = quiver!([xdose], [xLf-40], quiver=([0], [38]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([14500],[75],"16,600 g/m³")
    fig = quiver!([15000], [80], quiver=([1600], [35]), linecolor = :black, line=(:solid, 1))
    fig = plot!(
            xlabel = "Dose Concentration (g/m³)",
            ylabel = "Biofilm Thickness (μm)",
            #ylims = (0.0,200),
            xlims = (0,maximum(doses)),
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            legend = false,
            size = (800,500),
            margin = 10mm,
        )
    display(fig)
    savefig("Fig2a.svg")
end

begin # Plot Live vs doses (Average over biofilm)
    fig = plot()
    fig = plot!(doses,map(i -> meanLive(Pbs[i]),1:length(doses)),
        linewidth = 2,  
        linecolor = :black,    
    )
    #fig = annotate!([xdose],[15],@sprintf("%4.0f g/m³",xdose))
    #fig = quiver!([xdose], [20], quiver=([0], [16]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([14500],[5],"16,600 g/m³")
    fig = quiver!([15000], [10], quiver=([1600], [6]), linecolor = :black, line=(:solid, 1))
        
    fig = plot!(
        xlabel = "Dose Concentration (g/m³)",
        ylabel = "% Live",
        #ylims = (0.0,200),
        xlims = (0,maximum(doses)),
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=14,
        legend = false,
        size = (800,500),
        margin = 10mm,
    )
    display(fig)
    savefig("Fig2b.svg")
end

begin # Compute Log Reduction 
    for dose in [16_600] # dose(s) to compute log reduction at
        println(dose)
        # Find index closes to dose in doses 
        value,index_dose  =findmin(abs.(doses.-dose))
        value,index_nodose=findmin(abs.(doses.-0.0))
        # Compute # live cells
        Lf_nodose = Lfs[index_nodose][end]
        Lf_dose   = Lfs[index_dose  ][end]
        intVF_nodose = sum(Pbs[index_nodose][1,:])/length(Pbs[index_nodose][1,:])
        intVF_dose   = sum(Pbs[index_dose  ][1,:])/length(Pbs[index_dose  ][1,:])
        numlive_nodose = Lf_nodose*intVF_nodose
        numlive_dose   = Lf_dose  *intVF_dose  
        log_reduction = -log10(numlive_dose/numlive_nodose)
        println("Biofilm thickness")
        @printf("     no dose = %6.3f μm \n",Lf_nodose*1e6)
        @printf(" %6.0f dose = %6.3f μm \n", dose, Lf_dose*1e6)
        println("Mean volume fraction")
        @printf("     no dose = %6.3f \n",intVF_nodose)
        @printf(" %6.0f dose = %6.3f \n", dose, intVF_dose)
        @printf("Log reduction for %6.0f dose = -log₁₀( (%6.7f * %6.3f)/(%6.7f * %5.3f) \n"
            ,dose,Lf_dose,intVF_dose,Lf_nodose,intVF_nodose) 
        @printf("                              = %6.4f <=================\n",log_reduction)
    end
end


#########
# Fig 3: Concentration profiles within the biofilm. A, glucose concentration before and after biocide treatment. B, hydrogen peroxide concentration at steady state during continuous treatment.
#########
## Uses results from runs in Fig. 1 ##
begin # plot Glucose vs Depth in biofilm
    fig = plot()
    fig = plot!(zm1.*1e6,Sb1[1,:],line=(:solid, 2, :red),label="Dosing off")
    fig = plot!(zm2.*1e6,Sb2[1,:],line=(:dash,  2, :blue),label="Dosing on")
    fig = plot!(
        xlabel = "Height in Biofilm (μm)",
        ylabel = "Glucose Concentration (g/m³)",
        #ylims = (0.0,200),
        #xlims = (0,10),
        legend=:topleft,
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (800,500),
        margin = 10mm,
    )
    display(fig)
    savefig("Fig3a.svg")
end
begin # plot Hydrogen Peroxide vs Depth in biofilm
    fig = plot()
    fig = plot!(zm1.*1e6,Sb1[2,:],line=(:solid, 2, :red),label="Dosing off")
    fig = plot!(zm2.*1e6,Sb2[2,:],line=(:dash,  2, :blue),label="Dosing on")
    fig = plot!(
        xlabel = "Height in Biofilm (μm)",
        ylabel = "H₂O₂ Concentration (g/m³)",
        #ylims = (0.0,200),
        xlims = (0,maximum(zm2.*1e6)),
        legend=:topleft,
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (800,500),
        margin = 10mm,
    )
    display(fig)
    savefig("Fig3b.svg")
end

#########
# Fig 4: Spatial distribution of live and dead cells, A, and glucose consumption rate, B, within the biofilm before and after biocide treatment.
#########
## Uses results from runs in Fig. 1 ##
begin # plot live & dead vs Depth in biofilm
    fig = plot()
    fig = plot!(zm1.*1e6,Pb1[1,:],label="Live   - Dosing off",line=(:solid,2), linecolor = :red,)
    fig = plot!(zm1.*1e6,Pb1[2,:],label="Dead - Dosing off",  line=(:dash, 2), linecolor = :red,)
    fig = plot!(zm2.*1e6,Pb2[1,:],label="Live   - Dosing on ",line=(:solid,2), linecolor = :blue,)
    fig = plot!(zm2.*1e6,Pb2[2,:],label="Dead - Dosing on ",  line=(:dash, 2), linecolor = :blue,)
    fig = plot!(
        xlabel = "Height in Biofilm (μm)",
        ylabel = "Volume Fraction (-)",
        #ylims = (0.0,200),
        xlims = (0,maximum(zm2.*1e6)),
        legend=(0.91,0.95),
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (900,550),
        margin = 10mm,
        right_margin = 25mm,
    )
    display(fig)
    savefig("Fig4a.svg")
    
    println("Mean live volume fraction = $(meanLive(Pb2))")
    println("Mean Dead volume fraction = $(meanDead(Pb2))")
end
begin # plot Glucose Consumption Rate vs Depth in biofilm
    function consumption(Sb,Xb,Lf,t,zm,p)
        C = zeros(p.Nz)
        @unpack Nx,Nz,mu = p 
        μb=zeros(Nx,Nz)
        for j in 1:Nx
            C[:] = C[:] + mu[1](Sb,Xb,Lf,t,zm,p).*Xb[1,:]./p.Yxs[1,1] # Used by growth
        end
        return C
    end
    fig = plot()
    Xb1=similar(Pb1)
    Xb2=similar(Pb2)
    for j=1:p.Nx
        Xb1[j,:] = p.rho[j]*Pb1[j,:]
        Xb2[j,:] = p.rho[j]*Pb2[j,:]
    end
    fig = plot!(zm1.*1e6,consumption(Sb1,Pb1,Lf1,t1[end],zm1,p),line=(:solid, 2, :red), label="Dosing off")
    fig = plot!(zm2.*1e6,consumption(Sb2,Pb2,Lf2,t2[end],zm2,p),line=(:dash,  2, :blue),label="Dosing on")
    fig = plot!(
        xlabel = "Height in Biofilm (μm)",
        ylabel = "Glucose Consumption (g/m³-d)",
        #ylims = (0.0,200),
        xlims = (0,maximum(zm2.*1e6)),
        legend=:topleft,
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        foreground_color_legend = nothing,
        size = (800,500),
        margin = 10mm,
    )
    display(fig)
    savefig("Fig4b.svg")
end

#########
# Fig 5: Biofilm protection depends on hydrogen peroxide neutralizing activity of dead cells. A, biofilm thickness versus time with neutralization of hydrogen peroxide turned off for either live or dead cells. B, percent live cells versus time with neutralization of hydrogen peroxide turned off for either live or dead cells.
#This is relative to base case
#########
begin # run
    # Just run cases as sol is needed
    # if  isfile("cases4.jld2") && 
    #     isfile("cases5.jld2")
    #     println("Using saved files for cases4, and cases5")
    #     JLD2.@load "cases4.jld2" t4 zm4 X4 S4 Pb4 Sb4 Lf4 sol4
    #     JLD2.@load "cases5.jld2" t5 zm5 X5 S5 Pb5 Sb5 Lf5 sol5
    # else
        GlucoseIn = 100.0; # g/m³

        # Turn off live neutralization (only dead)
        dose1 = 500.0; dose2 = 0.0; k_bl=0.0; k_bd=10.0;
        t4,zm4,X4,S4,Pb4,Sb4,Lf4,sol4 = BiofilmSolver(p)

        # Turn off dead neutralization (only live)
        dose1 = 500.0; dose2 = 0.0; k_bl=10.0; k_bd=0.0;
        t5,zm5,X5,S5,Pb5,Sb5,Lf5,sol5 = BiofilmSolver(p)

    #    JLD2.@save "cases4.jld2" t4 zm4 X4 S4 Pb4 Sb4 Lf4 sol4
    #    JLD2.@save "cases5.jld2" t5 zm5 X5 S5 Pb5 Sb5 Lf5 sol5
    # end
end
begin # plot Biofilm Thickness vs Time
    fig = plot()
    fig = plot!(t2,Lf2'.*1e6,line=(:solid,2, :black), label="Live & Dead Neutralization")
    fig = plot!(t4,Lf4'.*1e6,line=(:dash, 2, :purple),label="Only Dead Neutralization")
    fig = plot!(t5,Lf5'.*1e6,line=(:dot,  2, :green), label="Only Live Neutralization")
    #fig = quiver!([2], [70], quiver=([0], [70]), linecolor = :black, line=(:solid, 1))
    #fig = annotate!([2],[60],"Dosing On")
    fig = quiver!([2], [20], quiver=([0], [40]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([1.9],[40],text("Dosing On",:right,12))
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
    display(fig)
    savefig("Fig5a.svg")
end
begin # plot % Live vs Time (average over biofilm)
    # Get mean substrates & particulates vs time 
    Sb_t2,Pb_t2=MeanBiofilmVarsWithTime(sol2,p)
    Sb_t4,Pb_t4=MeanBiofilmVarsWithTime(sol4,p)
    Sb_t5,Pb_t5=MeanBiofilmVarsWithTime(sol5,p)
    fig = plot()
    fig = plot!(t2,Pb_t2[1,:]./(Pb_t2[1,:].+Pb_t2[2,:]).*100,line=(:solid,2, :black), label="Live & Dead Neutralization")
    fig = plot!(t4,Pb_t4[1,:]./(Pb_t4[1,:].+Pb_t4[2,:]).*100,line=(:dash, 2, :purple),label="Only Dead Neutralization")
    fig = plot!(t5,Pb_t5[1,:]./(Pb_t5[1,:].+Pb_t5[2,:]).*100,line=(:dot,  2, :green), label="Only Live Neutralization")
    #fig = quiver!([2], [0], quiver=([0], [20]), linecolor = :black, line=(:solid, 1))
    #fig = annotate!([1.9],[10],text("Dosing On",:right,12))
    fig = quiver!([2], [10], quiver=([0], [20]), linecolor = :black, line=(:solid, 1))
    fig = annotate!([1.9],[20],text("Dosing On",:right,12))
    fig = plot!(
            xlabel = "Time (d)",
            ylabel = "% Live",
            #ylims = (0.0,200),
            xlims = (0,10),
            legend=(0.55,0.25),
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            foreground_color_legend = nothing,
            size = (800,500),
            margin = 10mm,
        )
    display(fig)
    savefig("Fig5b.svg")
end

#########
# Fig 6: Steady state biofilm thickness, A, and percent live cells, B, when subjected to continuous biocide treatment with inactivation of neutralizing capacity of either dead cells or live cells.
#### Same as Fig 2 with base-case (both live & dead , just live, just dead neutralization) ###
#########
begin # runs - no live neutralization
    if  isfile("case_noliveneut.jld2")
        println("Using saved files for case_noliveneut")
        JLD2.@load "case_noliveneut.jld2" ts_noliveneut zms_noliveneut Xs_noliveneut Ss_noliveneut Pbs_noliveneut Sbs_noliveneut Lfs_noliveneut
    else
        k_bl  =  0.0; #m³/g/d
        k_bd  = 10.0; #m³/g/d
        GlucoseIn = 100.0; # g/m³
        ts_noliveneut  = Vector{Float64}[]
        zms_noliveneut = Vector{Float64}[]
        Xs_noliveneut  = Matrix{Float64}[]
        Ss_noliveneut  = Matrix{Float64}[]
        Pbs_noliveneut = Matrix{Float64}[]
        Sbs_noliveneut = Matrix{Float64}[]
        Lfs_noliveneut = Matrix{Float64}[]
        for i in eachindex(doses)
            global dose1 = doses[i]
            global dose2 = 0.0
            println("Running solver with dose1=$dose1")
            t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p)
            push!( ts_noliveneut, t)
            push!(zms_noliveneut,collect(zm))
            push!( Xs_noliveneut, X)
            push!( Ss_noliveneut, S)
            push!(Pbs_noliveneut,Pb)
            push!(Sbs_noliveneut,Sb)
            push!(Lfs_noliveneut,Lf)
        end
        JLD2.@save "case_noliveneut.jld2" ts_noliveneut zms_noliveneut Xs_noliveneut Ss_noliveneut Pbs_noliveneut Sbs_noliveneut Lfs_noliveneut
    end
end
begin # runs - no dead neutralization
    if  isfile("case_nodeadneut.jld2")
        println("Using saved files for case_nodeadneut")
        JLD2.@load "case_nodeadneut.jld2" ts_nodeadneut zms_nodeadneut Xs_nodeadneut Ss_nodeadneut Pbs_nodeadneut Sbs_nodeadneut Lfs_nodeadneut
    else
        k_bl  = 10.0; #m³/g/d
        k_bd  =  0.0; #m³/g/d
        GlucoseIn = 100.0; # g/m³
        #doses = [1.0,10.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0,600.0,700.0,800.0,900.0,1000.0,2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,11000.0,12000.0,13000.0,14000.0,15000.0,16000.0,17000.0,18000.0,19000.0,20000.0]
        ts_nodeadneut  = Vector{Float64}[]
        zms_nodeadneut = Vector{Float64}[]
        Xs_nodeadneut  = Matrix{Float64}[]
        Ss_nodeadneut  = Matrix{Float64}[]
        Pbs_nodeadneut = Matrix{Float64}[]
        Sbs_nodeadneut = Matrix{Float64}[]
        Lfs_nodeadneut = Matrix{Float64}[]
        for i in eachindex(doses)
            global dose1 = doses[i]
            global dose2 = 0.0
            println("Running solver with dose1=$dose1")
            t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p)
            push!( ts_nodeadneut, t)
            push!(zms_nodeadneut,collect(zm))
            push!( Xs_nodeadneut, X)
            push!( Ss_nodeadneut, S)
            push!(Pbs_nodeadneut,Pb)
            push!(Sbs_nodeadneut,Sb)
            push!(Lfs_nodeadneut,Lf)
        end
        JLD2.@save "case_nodeadneut.jld2" ts_nodeadneut zms_nodeadneut Xs_nodeadneut Ss_nodeadneut Pbs_nodeadneut Sbs_nodeadneut Lfs_nodeadneut
    end
end

begin # Plot Lf vs doses 
    fig = plot()
    fig = plot!(doses,map(i -> Lfs[i][end]*1e6,1:length(doses)),           line=(:solid, 2, :black), label="Live & Dead Neutralization")
    fig = plot!(doses,map(i -> Lfs_noliveneut[i][end]*1e6,1:length(doses)),line=(:dot,   2, :purple),label="Only Dead Neutralization")
    fig = plot!(doses,map(i -> Lfs_nodeadneut[i][end]*1e6,1:length(doses)),line=(:dash,  2, :green), label="Only Live Neutralization")
    #fig = plot!([0,maximum(doses)],[Lf1[end],Lf1[end]].*1e6)
    fig = plot!(
            xlabel = "Dose (g/m³)",
            ylabel = "Biofilm Thickness (μm)",
            #ylims = (0.0,200),
            xlims = (0,maximum(doses)),
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            legend = :topright,
            foreground_color_legend = nothing,  
            size = (800,500),
            margin = 10mm,
        )
    display(fig)
    savefig("Fig6a.svg")
    for i=1:5
        @printf("%6.0f  %6.3f \n",doses[i],Lfs_nodeadneut[i][end]*1e6)
    end
end
begin # Plot % Live vs doses
    fig = plot()
    fig = plot!(doses,map(i -> meanLive(Pbs[i]),1:length(doses))           ,line=(:solid, 2, :black), label="Live & Dead Neutralization")
    fig = plot!(doses,map(i -> meanLive(Pbs_noliveneut[i]),1:length(doses)),line=(:dot,   2, :purple),label="Only Dead Neutralization")
    fig = plot!(doses,map(i -> meanLive(Pbs_nodeadneut[i]),1:length(doses)),line=(:dash,  2, :green), label="Only Live Neutralization")
    #fig = plot!([0,maximum(doses)],[Lf1[end],Lf1[end]].*1e6)
    fig = plot!(
            xlabel = "Dose (g/m³)",
            ylabel = "% Live",
            #ylims = (0.0,200),
            xlims = (0,maximum(doses)),
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            legend = :topright,
            foreground_color_legend = nothing,  
            size = (800,500),
            margin = 10mm,
        )
    display(fig)
    savefig("Fig6b.svg")
end

#########
# Fig 7: Biofilm tolerance of HP depends on glucose concentration. A, steady state biofilm thickness when subjected to continuous HP treatment at various glucose concentrations. B, steady state percent live cells when subjected to continuous HP treatment at various glucose concentrations.
#########
begin # runs (include lower concentrations - expect thickness to drop to zero at some low enough conc.)
    #GlucoseIns = [5,10,20,30,40,50.0,60,70,80,90,100.0,120,130,140,150.0,160,170,180,190,200.0]; # g/m³
    GlucoseIns = 0.0:1.0:200.0
    if  isfile("case_g.jld2")
        println("Using saved files for case_g")
        JLD2.@load "case_g.jld2" t_g zm_g X_g S_g Pb_g Sb_g Lf_g
    else
        k_bl  = 10.0; #m³/g/d
        k_bd  = 10.0; #m³/g/d
        dose1 = 500; # g/m³
        dose2 = 0.0; # g/m³
        t_g  = Vector{Float64}[]
        zm_g = Vector{Float64}[]
        X_g  = Matrix{Float64}[]
        S_g  = Matrix{Float64}[]
        Pb_g = Matrix{Float64}[]
        Sb_g = Matrix{Float64}[]
        Lf_g = Matrix{Float64}[]
        for i in eachindex(GlucoseIns)
            global GlucoseIn = GlucoseIns[i]
            t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p)
            push!( t_g, t)
            push!(zm_g,collect(zm))
            push!( X_g, X)
            push!( S_g, S)
            push!(Pb_g,Pb)
            push!(Sb_g,Sb)
            push!(Lf_g,Lf)
        end
        JLD2.@save "case_g.jld2" t_g zm_g X_g S_g Pb_g Sb_g Lf_g
    end
end

begin # Plot Lf vs GlucoseIn 
    fig = plot()
    fig = plot!(GlucoseIns,map(i -> Lf_g[i][end]*1e6,1:length(GlucoseIns)),line=(:solid, 2, :black))
    fig = plot!(
            xlabel = "Glucose Influent Concentration (g/m³)",
            ylabel = "Biofilm Thickness (μm)",
            #ylims = (0.0,200),
            xlims = (0,maximum(GlucoseIns)),
            xguidefontsize=16,
            yguidefontsize=16,
            xtickfontsize=14,
            ytickfontsize=14,
            legendfontsize=12,
            legend = false,
            size = (800,500),
            margin = 10mm,
        )
    display(fig)
    savefig("Fig7a.svg")
    for i=1:20
        @printf(" %6.3f  %6.3f \n", GlucoseIns[i],Lf_g[i][end]*1e6)
    end
end
begin # Plot Live vs GlucoseIn (need depth averaged)
    plot()
    plot!(GlucoseIns,map(i -> meanLive(Pb_g[i]),1:length(GlucoseIns)),line=(:solid, 2, :black))
    fig = plot!(
        xlabel = "Glucose Influent Concentration (g/m³)",
        ylabel = "% Live",
        #ylims = (0.0,200),
        xlims = (0,maximum(GlucoseIns)),
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        legend = false,
        size = (800,500),
        margin = 10mm,
    )
    display(fig)
    savefig("Fig7b.svg")
end

