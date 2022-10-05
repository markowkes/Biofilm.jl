using Plots
using LaTeXStrings
using Measures
using ProgressMeter
using UnPack
using Printf

# Called during ODE solve to see solution progress
function outputs(integrator)

    # Unpack integrator
    sol=integrator.sol 
    p=integrator.p[1]
    r=integrator.p[2]
    @unpack Nz,outPeriod,plotPeriod= p

    # Perform plots on output period 
    modt=mod(sol.t[end],outPeriod)
    if modt≈0.0 || modt≈outPeriod


        # Convert solution to dependent variables
        t,X,S,Pb,Sb,Lf=unpack_solutionForPlot(sol,p,r)

        # Print titles to REPL every 10 outPeriod
        if mod(sol.t[end],outPeriod*10)≈0.0 || mod(sol.t[end],outPeriod*10)≈outPeriod*10
            printBiofilmTitles(p)
        end

        # Print values to REPL every 1 outPeriod
        printBiofilmValues(t[end],X[:,end],S[:,end],Pb,Sb,Lf[end],p)

        # Plot results
        if p.makePlots && (mod(sol.t[end],plotPeriod)≈0.0 || mod(sol.t[end],plotPeriod)≈plotPeriod)
            makePlots(t,X,S,Pb,Sb,Lf,p)
        end

    end
end 

# Take in array and compute reasonable ylimits
function pad_ylim(A)
    ymin=minimum(A)
    ymax=maximum(A)
    yavg=0.5*(ymin+ymax)
    deltay = max(0.1*yavg,0.6*(ymax-ymin))
    return [max(0.0,yavg-deltay),yavg+deltay]
end 

# Make plots
function makePlots(t,X,S,Pb,Sb,Lf,p)
    @unpack Nx,Ns,Nz,Title,XNames,SNames,Ptot,rho = p 

    # Adjust names to work with legends
    Nx==1 ? Xs=XNames[1] : Xs=reshape(XNames,1,length(XNames))
    Ns==1 ? Ss=SNames[1] : Ss=reshape(SNames,1,length(SNames))

    # Compute grid
    z=range(0.0,Lf[end],Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])
    dz=z[2]-z[1]
    g=biofilmGrid(z,zm,dz)

    # Tank particulate concentration
    p1=plot(t,X',label=Xs,ylim=pad_ylim(X))
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Tank~Biomass~Concentration~} [g/m^3]")

    # Tank substrate concentration
    p2=plot(t,S',label=Ss,ylim=pad_ylim(S))
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Tank~Substrate~Concentration~} [g/m^3]")

    # Biofilm thickness
    p3=plot(t,Lf'*1e6,label="Thickness",ylim=pad_ylim(Lf*1e6))
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Biofilm~Thickness~} [μm]")

    # Biofilm particulate volume fractioin 
    p4=plot(zm,Pb',label=Xs,ylim=pad_ylim(Pb))
    #p4=plot!(zm,sum(Pb,dims=1)',label="Sum")
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Biofilm~Particulate~Volume~Fraction~[-]}")
    
    # Biofilm substrate concentration
    p5=plot(zm,Sb',label=Ss,ylim=pad_ylim(Sb))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Biofilm~Substrate~Concentration~} [g/m^3]")

    # Particulate growthrates vs depth
    Xb=similar(Pb)
    for j=1:Nx
        Xb[j,:] = rho[j]*Pb[j,:]  # Compute particulate concentrations
    end
    μb    = computeMu_biofilm(Sb,Xb,Lf[end],t[end],p,g)   # Growthrates in biofilm
    p6=plot(zm,μb',label=Xs,ylim=pad_ylim(μb))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Particulate~Growthrates~[-]}")

    # Put plots together
    myplt=plot(p1,p2,p3,p4,p5,p6,
        layout=(2,3),
        size=(1600,1000),
        plot_title=@sprintf("%s : t = %.2f",Title,t[end]),
        #plot_titlevspan=0.5,
        left_margin=10mm, 
        bottom_margin=10mm,
        foreground_color_legend = nothing,
        legend = :outertop,
    )
    display(myplt)
    return
end

# Make plots
function makeBiofilmPlots(t,Pb,Sb,Lf,p)
    @unpack Nx,Ns,Nz,Title,XNames,SNames,Ptot,mu,rho = p 

    # Adjust names to work with legends
    Nx==1 ? Xs=XNames[1] : Xs=reshape(XNames,1,length(XNames))
    Ns==1 ? Ss=SNames[1] : Ss=reshape(SNames,1,length(SNames))

    # Package grid
    z=range(0.0,Lf[end],Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])
    dz=z[2]-z[1]
    g=biofilmGrid(z,zm,dz)

    # Make plots
    p1=plot(zm,Pb',label=Xs,ylim=pad_ylim(Pb))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Particulate~Volume~Fraction~[-]}")
    
    p2=plot(zm,Sb',label=Ss,ylim=pad_ylim(Sb))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Substrate~Concentration~} [g/m^3]")

    # Particulate growthrates vs depth
    Xb=similar(Pb)
    for j=1:Nx
        Xb[j,:] = rho[j]*Pb[j,:]  # Compute particulate concentrations
    end
    μb    = computeMu_biofilm(Sb,Xb,Lf,t,p,g)   # Growthrates in biofilm
    p3=plot(zm,μb',label=Xs,ylim=pad_ylim(μb))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Particulate~Growthrates~[-]}")
    

    # Put plots together
    myplt=plot(p1,p2,p3,
        layout=(1,3),
        size=(1600,500),
        plot_title=@sprintf("%s : t = %.2f",Title,t[end]),
        #plot_titlevspan=0.5,
        left_margin=10mm, 
        bottom_margin=10mm,
        foreground_color_legend = nothing,
    )
    display(myplt)
    return
end