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

    # Perform plots on output period 
    modt=mod(sol.t[end],p.outPeriod)
    if modt≈0.0 || modt≈p.outPeriod

        # Convert solution to dependent variables
        t,X,S,Pb,Sb,Lf=unpack_solution(sol,p,r)

        println("Time = ",sol.t[end])   

        # Plot results
        if p.makePlots 
            outputs(t,X,S,Pb,Sb,Lf,p)
        end

    end
end 

# Take in array and compute reasonable ylimits
function pad_ylim(A)
    ymin=minimum(A)
    ymax=maximum(A)
    yavg=0.5*(ymin+ymax)
    deltay = max(0.1*yavg,0.5*(ymax-ymin))
    return [yavg-deltay,yavg+deltay]
end 

# Make plots
function outputs(t,X,S,Pb,Sb,Lf,p)
    @unpack Nx,Ns,Nz,Title,XNames,SNames,Ptot = p 

    # Adjust names to work with legends
    Nx==1 ? Xs=XNames[1] : Xs=reshape(XNames,1,length(XNames))
    Ns==1 ? Ss=SNames[1] : Ss=reshape(SNames,1,length(SNames))

    # Final biofilm grid
    z=range(0.0,Lf[end],Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])

    # Make plots
    p1=plot(t,X',label=Xs,ylim=pad_ylim(X))
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Biomass~Concentration~} [g/m^3]")

    p2=plot(t,S',label=Ss,ylim=pad_ylim(S))
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Substrate~Concentration~} [g/m^3]")

    p3=plot(t,Lf'*1e6,legend=false,ylim=pad_ylim(Lf*1e6))
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Thickness~} [μm]")

    p4=plot(zm,Pb',label=Xs,ylim=(0,min(1,Ptot+0.1)))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Particulate~Volume~Fraction~[-]}")
    
    p5=plot(zm,Sb',label=Ss,ylim=pad_ylim(Sb))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Substrate~Concentration~} [g/m^3]")

    # Put plots together
    myplt=plot(p1,p2,p3,p4,p5,
        layout=(2,3),
        size=(1600,1000),
        plot_title=@sprintf("%s : t = %.2f",Title,t[end]),
        #plot_titlevspan=0.5,
        left_margin=10mm, 
        bottom_margin=10mm,
    )
    display(myplt)
    return
end