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

    # Convert solution to dependent variables
    t,X,S,Pb,Sb,Lf=unpack_solution(sol,p,r)

    println("Time = ",integrator.sol.t[end])   

    # Plot results
    if p.makePlots 
        outputs(t,X,S,Pb,Sb,Lf,p)
    end
end 

function outputs(t,X,S,Pb,Sb,Lf,p)
    @unpack Nx,Ns,Nz,Title,XNames,SNames,Ptot = p 

    # Adjust names to work with legends
    Nx==1 ? Xs=XNames[1] : Xs=reshape(XNames,1,length(XNames))
    Ns==1 ? Ss=SNames[1] : Ss=reshape(SNames,1,length(SNames))

    # Final biofilm grid
    z=range(0.0,Lf[end],Nz+1)
    zm=0.5*(z[1:Nz]+z[2:Nz+1])

    # Make plots
    p1=plot(t,X',label=Xs)
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Biomass Con.~} [g/m^3]")

    p2=plot(t,S',label=Ss)
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Substrate Con.~} [g/m^3]")

    p3=plot(t,Lf'*1e6,legend=false)
    xaxis!(L"\textrm{Time~[days]}")
    yaxis!(L"\textrm{Thickness~} [μm]")

    p4=plot(zm,Pb',label=Xs,ylim=(0,min(1,Ptot+0.1)))
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Particulate~Volume~Fraction~[-]}")
    
    p5=plot(zm,Sb',label=Ss)
    xaxis!(L"\textrm{Thickness~} [\mu m]")
    yaxis!(L"\textrm{Substrate~Concentration~} [g/m^3]")

    # N=length(t)
    # dt=t[2:N]-t[1:N-1]
    # #display(dt)
    # p6=plot(1:(N-1),dt)
    # xaxis!(L"\textrm{Iteration }")
    # yaxis!(L"\textrm{Timestep } [s]")


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