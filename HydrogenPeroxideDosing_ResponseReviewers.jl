

begin # Run HydrogenPeroxideDosing.jl 
    include("HydrogenPeroxideDosing.jl")
end

# ---------------------
# Reviewer 1, Comment 1
# ---------------------
begin # Vary k_bd from 0 to 10 m^3/g/d and compute L_f and % Live
    k_bds = range(0,10,20)
    thickness=zeros(size(k_bds))
    perLive  =zeros(size(k_bds))
    for n in 1:length(k_bds)
        k_bl  = 10.0; #m³/g/d
        k_bd  = k_bds[n]
        GlucoseIn = 100.0; # g/m³
        println("Running with k_bd=$k_bd")
        dose1 = 500.0; dose2 =    0.0; t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p)
        thickness[n] = Lf[end]
        Sb_t,Pb_t=MeanBiofilmVarsWithTime(sol,p)
        perLive[n]   = Pb_t[1,end]./(Pb_t[1,end].+Pb_t[2,end]).*100
    end
end
begin # Plot
    fig = plot(k_bds,thickness*1e6,label="")
    fig = plot!(
        xlabel = "Biocide Neutralization Rate of Dead Biomass",
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
    savefig("reviewer1_comment1.pdf")
end

# ---------------------
# Reviewer 2 - Comment on need for L_L 
# ---------------------
begin # Run with L_L = 0 and compare with run 2 
    k_bl  = 10.0; #m³/g/d
    k_bd  = 10.0; #m³/g/d
    GlucoseIn = 100.0; # g/m³
    p_noLL = @set p.LL = 0.0
    dose1 = 500.0; dose2 =    0.0; t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p_noLL)
end
begin # Plot
    fig = plot()
    fig = plot!(t2,1e6.*Lf2',label="L_L = $(p.LL)")
    fig = plot!(t ,1e6.*Lf' ,label="L_L = $(p_noLL.LL)")
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
    display(fig)
    savefig("reviewer2_noLL.pdf")
end

# ---------------------
# Reviewer 2 - Longer run times to get Fig 4 to steady-state 
# ---------------------
begin # Runs with varying tFinal
    for tFinal in [] #[] #[10,100,1000,2000,5000,10000] # takes a while to run - so comment out to reproduce
        println("Working on tFinal = $tFinal")
        k_bl  = 10.0; #m³/g/d
        k_bd  = 10.0; #m³/g/d
        GlucoseIn = 100.0; # g/m³
        p_tFinal = @set p.tFinal = tFinal
        dose1 =   0.0; dose2 =    0.0; t1vt,zm1vt,X1vt,S1vt,Pb1vt,Sb1vt,Lf1vt,sol1vt = BiofilmSolver(p_tFinal)
        dose1 = 500.0; dose2 =    0.0; t2vt,zm2vt,X2vt,S2vt,Pb2vt,Sb2vt,Lf2vt,sol2vt = BiofilmSolver(p_tFinal)

        # Plot volume fraction
        fig = plot()
        fig = plot!(zm1vt.*1e6,Pb1vt[1,:],label="Live   - No Dosing",line=(:solid,2), linecolor = :red,)
        fig = plot!(zm1vt.*1e6,Pb1vt[2,:],label="Dead - No Dosing",  line=(:dash, 2), linecolor = :red,)
        fig = plot!(zm2vt.*1e6,Pb2vt[1,:],label="Live   - Dosing On ",line=(:solid,2), linecolor = :blue,)
        fig = plot!(zm2vt.*1e6,Pb2vt[2,:],label="Dead - Dosing On ",  line=(:dash, 2), linecolor = :blue,)
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
        savefig("reviewer2_steadyState_$(tFinal).pdf")
    end
end

begin # Runs with tFinal = 2000 to make results for paper
    tFinal = 2000
    k_bl  = 10.0; #m³/g/d
    k_bd  = 10.0; #m³/g/d
    GlucoseIn = 100.0; # g/m³
    p_tFinal = @set p.tFinal = tFinal
    dose1 =   0.0; dose2 =    0.0; t1vt,zm1vt,X1vt,S1vt,Pb1vt,Sb1vt,Lf1vt,sol1vt = BiofilmSolver(p_tFinal)
    dose1 = 500.0; dose2 =    0.0; t2vt,zm2vt,X2vt,S2vt,Pb2vt,Sb2vt,Lf2vt,sol2vt = BiofilmSolver(p_tFinal)

    # Plot volume fraction
    fig = plot()
    fig = plot!(zm1vt.*1e6,Pb1vt[1,:],label="Live   - No Dosing",line=(:solid,2), linecolor = :red,)
    fig = plot!(zm1vt.*1e6,Pb1vt[2,:],label="Dead - No Dosing",  line=(:dash, 2), linecolor = :red,)
    fig = plot!(zm2vt.*1e6,Pb2vt[1,:],label="Live   - Dosing On ",line=(:solid,2), linecolor = :blue,)
    fig = plot!(zm2vt.*1e6,Pb2vt[2,:],label="Dead - Dosing On ",  line=(:dash, 2), linecolor = :blue,)
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
    savefig("Fig4a_$(tFinal)days.pdf")
    println("Mean live volume fraction = $(meanLive(Pb2vt))")
    println("Mean Dead volume fraction = $(meanDead(Pb2vt))")
end

# ---------------------
# Reviewer 3 - Sensitivity of 16,600 g/m^3 dosing
# ---------------------
# Vary initial conditions to see if Fig 2a changes (focusing on area near cutoff)
begin # runs
    Lfends = []
    Lfos=[5e-6, 50.0e-6, 500.0e-6]
    for Lfo in Lfos
    p_thickness = @set p.Lfo = Lfo
    doses = 16000.0:200.0:17400.0
    # if  isfile("doses_inputs.jld2") # Check for file of saved results to reduce runtime
    #     println("Using saved files for doses")
    #     JLD2.@load "doses.jld2" ts zms Xs Ss Pbs Sbs Lfs
    # else # Run simulations and save results 
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
        Lfend = []
        for i in eachindex(doses)
            global dose1 = doses[i]
            global dose2 = 0.0
            println("Running solver with dose1=$dose1")
            t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p_thickness)
            push!( ts, t)
            push!(zms,collect(zm))
            push!( Xs, X)
            push!( Ss, S)
            push!(Pbs,Pb)
            push!(Sbs,Sb)
            push!(Lfs,Lf)
            push!(Lfend,Lf[end])
        end
        push!(Lfends,Lfend)
        # JLD2.@save "doses.jld2" ts zms Xs Ss Pbs Sbs Lfs
    
        # Plot Lf vs doses (add additional doses - make sure at S.S.)
        fig = plot()
        # fig = plot!([0,maximum(doses)],[Lf1[end],Lf1[end]].*1e6, 
        #     linewidth = 1, 
        #     linecolor = :black,
        #     line=(:dot),
        # )
        #fig = annotate!(1.5e4,150,"No Dose Thickness")
        fig = plot!(doses,map(i -> Lfs[i][end]*1e6,1:length(doses)),
            linewidth = 2,
            #linecolor = :black,
        )
        # Add label at crossover point 
        # A = map(i -> Lfs[i][end]*1e6,1:length(doses))
        # A = abs.(A .- A[1])
        # A[1]=100.0
        # value,index = findmin(A)
        # xdose = doses[index]
        # xLf   = Lfs[index][end]*1e6
        # fig = annotate!([xdose],[xLf-50],@sprintf("%4.0f g/m³",xdose))
        # fig = quiver!([xdose], [xLf-40], quiver=([0], [38]), linecolor = :black, line=(:solid, 1))
        # fig = annotate!([14500],[75],"16,600 g/m³")
        # fig = quiver!([15000], [80], quiver=([1600], [35]), linecolor = :black, line=(:solid, 1))
        fig = plot!(
                xlabel = "Dose Concentration (g/m³)",
                ylabel = "Biofilm Thickness (μm)",
                #ylims = (0.0,200),
                xlims = (minimum(doses),maximum(doses)),
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
        savefig("Fig2a_Lfo$Lfo.pdf")
    end
end
begin # Plot and analysis
    fig = plot()
    for i=1:length(Lfos)
        fig = plot!( doses,1e6.*Lfends[i],label="$(Lfos[i])")
        # Percent diff from Lfo=50e-6 (baseline case)
        perDiff = mean( (Lfends[i].-Lfends[2])./Lfends[2].*100.0 )
        println("$(Lfos[i]) mean(% diff) = $perDiff")
    end
    fig = plot!(
        xlabel = "Dose Concentration (g/m³)",
        ylabel = "Biofilm Thickness (μm)",
        xlims = (minimum(doses),maximum(doses)),
        xguidefontsize=16,
        yguidefontsize=16,
        xtickfontsize=14,
        ytickfontsize=14,
        legendfontsize=12,
        #legend = false,
        size = (800,500),
        margin = 10mm,
    )
    display(fig)
    savefig("reviewer3_varyingIC.pdf")

end