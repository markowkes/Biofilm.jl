using Plots
using LaTeXStrings
using Measures
using ProgressMeter
using UnPack
using Printf

""" 
    outputs(integrator)
Called during ODE solve to see solution progress
"""
function outputs(integrator)

    # Unpack integrator
    sol=integrator.sol 
    p=integrator.p
    @unpack Nz,outPeriod,plotPeriod= p

    # Perform plots on output period 
    modt=mod(sol.t[end],outPeriod)
    if modt≈0.0 || modt≈outPeriod


        # Convert solution to dependent variables
        t,Xt,St,Pb,Sb,Lf=unpack_solutionForPlot(sol,p)

        # Print titles to REPL every 10 outPeriod
        if mod(sol.t[end],outPeriod*10)≈0.0 || mod(sol.t[end],outPeriod*10)≈outPeriod*10
            printBiofilmTitles(p)
        end

        # Print values to REPL every 1 outPeriod
        printBiofilmValues(t[end],Xt[:,end],St[:,end],Pb,Sb,Lf[end],p)

        # Plot results
        if p.makePlots && (mod(sol.t[end],plotPeriod)≈0.0 || mod(sol.t[end],plotPeriod)≈plotPeriod)
            plt = biofilm_plot(sol,p)
            display(plt)
        end

    end
end 

"""
    pad_ylim(A)
Take in array A and compute reasonable ylimits ylmin = pad_ylmin(A)
""" 
function pad_ylim(A)
    ymin=minimum(A)
    ymax=maximum(A)
    yavg=0.5*(ymin+ymax)
    deltay = max(0.1*yavg,0.6*(ymax-ymin))
    return [max(0.0,yavg-deltay),yavg+deltay]
end 

""" 
    processRecipeInputs(h)
Process inputs to biofilm_plot plot recipe 
"""
function processRecipeInputs(h) 

    # Check first argument
    if length(h.args) >= 2
        # Get & check solution 
        sol = h.args[1]
        try 
            sol.t 
            sol.u
        catch
            error("1st argument to biofilm_plot should be solution, e.g., biofilm_plot(sol)")
        end

        # Get param 
        p = h.args[2]
        has_p = true
    else
        error("biofilm_plot requries at least two inputs, the solution and parameters e.g., biofilm_plot(sol,p)") 
    end
    
    # Check if third argument is present 
    if length(h.args) >= 3
        if typeof(h.args[3]) <: String
            desc = " : " * h.args[3]
        else
            error("3nd argument to biofilm_plot should be string.  Got: $(typeof(h.args[3]))")
        end
    else
        desc = ""
    end

    # Check paramters and add defaults 
    p = checkInputs(p)

    # Unpack some parameters
    @unpack Nx,Ns,XNames, SNames = p

    # Compute ranges of dependent variables in solution vector 
    p = computeRanges(p)

    # Extract variables from solution 
    t,Xt,St,Pb,Sb,Lf = unpack_solutionForPlot(sol,p)

    return p,desc,t,Xt,St,Pb,Sb,Lf
end

"""
    biofilm_plot(sol,p)
    biofilm_plot(sol,p,desc)

Plot the solution of a biofilm simulation.  Makes 6 plots Xt(t), St(t), Lf(t), Pb(z), Sb(z), μ(z) or src(z).  
"""
@userplot Biofilm_Plot
@recipe function f(h::Biofilm_Plot)

    p,desc,t,Xt,St,Pb,Sb,Lf = processRecipeInputs(h)

    # set up the subplots 
    layout --> (2,3)
    left_margin   --> 10mm 
    bottom_margin --> 10mm
    foreground_color_legend --> nothing
    legend --> :outertop
    size --> p.plotSize
    plot_title := @sprintf("%s : t = %.2f",p.Title,t[end])

    # Create grid for plots within biofilm 
    g = computeGrid(Lf[end],p)
    zm = g.zm
    
    # Tank particulate concentration
    for n in eachindex(view(Xt,:,1))
        @series begin 
            subplot := 1
            label := p.XNames[n] .* desc
            xaxis := (L"\textrm{Time~[days]}")
            yaxis := (L"\textrm{Tank~Particulate~Conc.~} [g/m^3]")
            x=t 
            y=Xt[n,:]
            ylim := pad_ylim(Xt)
            x,y
        end
    end    
    
    # Tank solute concentration
    for n in eachindex(view(St,:,1))
        @series begin 
            subplot := 2
            label := p.SNames[n] .* desc
            xaxis := (L"\textrm{Time~[days]}")
            yaxis := (L"\textrm{Tank~Solute~Conc.~} [g/m^3]")
            x=t 
            y=St[n,:]
            ylim := pad_ylim(St)
            (x,y)
        end
    end

    # Biofilm thickness
    @series begin
        subplot := 3
        label := "Thickness" .* desc
        xaxis := (L"\textrm{Time~[days]}")
        yaxis := (L"\textrm{Biofilm~Thickness~} [μm]")
        x = t
        y = Lf'.*1e6
        ylim := pad_ylim(y)
        x,y
    end
    
    # Biofilm particulate volume fraction 
    for n in eachindex(view(Pb,:,1))
        @series begin
            subplot := 4
            label := p.XNames[n] .* desc
            xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
            yaxis := (L"\textrm{Biofilm~Particulate~Vol.~Frac.~[-]}")
            x = 1e6.*zm
            y = Pb[n,:]
            ylim := pad_ylim(Pb)
            x, y
        end
    end

    # Biofilm solute concentration
    for n in eachindex(view(Sb,:,1))
        @series begin
            subplot := 5
            label := p.SNames[n] .* desc
            xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
            yaxis := (L"\textrm{Biofilm~Solute~Conc.~} [g/m^3]")   
            x = 1e6.*zm
            y = Sb[n,:]
            ylim := pad_ylim(Sb)
            x, y
        end
    end

    # Optional 6th plot
    @unpack Nx,Nz,optionalPlot,rho,srcX = p
    if optionalPlot == "growthrate"
        # Particulate growthrates vs depth
        Xb=similar(Pb)
        for j=1:Nx
            Xb[j,:] = rho[j]*Pb[j,:]  # Compute particulate concentrations
        end
        μb    = computeMu_biofilm(Sb,Xb,Lf[end],t[end],p,g)   # Growthrates in biofilm
        for n in eachindex(view(μb,:,1))
            @series begin 
                subplot := 6    
                xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
                yaxis := (L"\textrm{Particulate~Growthrates~[-]}")
                x = 1e6.*zm 
                y = μb[n,:]
                label := p.XNames[n] .* desc
                ylim := pad_ylim(μb)
                x,y 
            end
        end
        
    elseif optionalPlot == "source"
        # Particulate source term vs depth
        srcs=similar(Pb)
        for i=1:Nz
            for j=1:Nx
                srcs[j,i]=srcX[j](Sb[:,i],Pb[:,i]*rho[j],Lf[end],t,zm[i],p)[1]
            end
        end
        for n in eachindex(view(srcs,:,1))
            @series begin
                subplot := 6
                xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
                yaxis := (L"\textrm{Particulate~Source~} [g/m^3\cdot d]")
                x = 1e6.*zm
                y = srcs[n,:]
                label := p.XNames[n] .* desc
                x,y
            end
        end
    else
        error("Unknown optionalPlot type of $optionalPlot.  Should be \"growthrate\" or \"source\".")
    end
end


"""
    biofilm_plot_Film(sol,p)
    biofilm_plot_Film(sol,p,desc)

Plot the solution of a biofilm simulation
"""
@userplot Biofilm_Plot_Film
@recipe function f(h::Biofilm_Plot_Film)

    p,desc,t,Xt,St,Pb,Sb,Lf = processRecipeInputs(h)

    # set up the subplots 
    layout --> (1,3)
    left_margin   --> 10mm 
    bottom_margin --> 10mm
    foreground_color_legend --> nothing
    legend --> :outertop
    size --> p.plotSize
    plot_title := @sprintf("%s : t = %.2f",p.Title,t[end])

    # Create grid for plots within biofilm 
    g = computeGrid(Lf[end],p)
    zm = g.zm

    # Biofilm particulate volume fraction 
    for n in eachindex(view(Pb,:,1))
        @series begin
            subplot := 1
            label := p.XNames[n] .* desc
            xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
            yaxis := (L"\textrm{Biofilm~Particulate~Vol.~Frac.~[-]}")
            x = 1e6.*zm
            y = Pb[n,:]
            ylim := pad_ylim(Pb)
            x, y
        end
    end

    # Biofilm solute concentration
    for n in eachindex(view(Sb,:,1))
        @series begin
            subplot := 2
            label := p.SNames[n] .* desc
            xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
            yaxis := (L"\textrm{Biofilm~Solute~Conc.~} [g/m^3]")   
            x = 1e6.*zm
            y = Sb[n,:]
            ylim := pad_ylim(Sb)
            x, y
        end
    end

    # Optional 6th plot
    @unpack Nx,Nz,optionalPlot,rho,srcX = p
    if optionalPlot == "growthrate"
        # Particulate growthrates vs depth
        Xb=similar(Pb)
        for j=1:Nx
            Xb[j,:] = rho[j]*Pb[j,:]  # Compute particulate concentrations
        end
        μb    = computeMu_biofilm(Sb,Xb,Lf[end],t[end],p,g)   # Growthrates in biofilm
        for n in eachindex(view(μb,:,1))
            @series begin 
                subplot := 3
                xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
                yaxis := (L"\textrm{Particulate~Growthrates~[-]}")
                x = 1e6.*zm 
                y = μb[n,:]
                label := p.XNames[n] .* desc
                ylim := pad_ylim(μb)
                x,y 
            end
        end
        
    elseif optionalPlot == "source"
        # Particulate source term vs depth
        srcs=similar(Pb)
        for i=1:Nz
            for j=1:Nx
                srcs[j,i]=srcX[j](Sb[:,i],Pb[:,i]*rho[j],Lf[end],t,zm[i],p)[1]
            end
        end
        for n in eachindex(view(srcs,:,1))
            @series begin
                subplot := 3
                xaxis := (L"\textrm{Height~in~Biofilm~} [\mu m]")
                yaxis := (L"\textrm{Particulate~Source~} [g/m^3\cdot d]")
                x = 1e6.*zm
                y = srcs[n,:]
                label := p.XNames[n] .* desc
                x,y
            end
        end
    else
        error("Unknown optionalPlot type of $optionalPlot.  Should be \"growthrate\" or \"source\".")
    end
end