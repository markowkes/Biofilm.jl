# test_Diffusion.jl
# -----------------------------------------------
# Checks the rate of diffusion of a solute
# into the biofilm.  Runs simulations on different 
# grids and checks the convergence rate 
# -----------------------------------------------

using Biofilm
using Statistics
using Printf
# using Plots

function test_Diffusion()

    # Loop over # grid points
    Nzs = [10,25,50,75,100]
    error = similar(Nzs,Float64)
    for n in eachindex(Nzs)

        # Constants used for growthrates of particulate(s)
        mumax = 2000;
        KM = 2500;

        # Define a structure to hold all the parameters
        p = (
            # --------------------- #
            # Simulation Parameters #
            # --------------------- #
            Title="Diffusion Test",
            tFinal=1,       # Simulation time [days]
            tol=1e-4,       # Tolerance
            outPeriod=1.0,  # Time between outputs [days]
            makePlots=false,

            # ---------------------- #
            # Particulate Parameters #
            # ---------------------- #
            XNames=["Bug"],     # Particulate names
            Xto=[10.0],          # Tank particulate concentraion initial condition(s)
            Pbo=[0.08],         # Biofilm particulates volume fraction initial condition(s) 
            rho=[2.0E4],        # Particulate densities
            Kdet=20000.0,       # Particulates detachment coefficient
            srcX=[(S,X,Lf,t,z,p) -> 0.0],     # Source of particulates
            # Growthrates for each particulate (constants defined above!)
            mu=[(S,X,Lf,t,z,p) -> (mumax * S[1]) ./ (KM)],

            # -------------------- #
            # Solute Parameters    #
            # -------------------- #
            SNames=["Oxygen"],   # Solute names
            Sin=[(t) -> 100],    # Solute inflow (can be function of time)
            Sto=[25.0],          # Tank solute concentraion initial condition(s)
            Sbo=[0.0],           # Biofilm solutes concentration initial condition(s)
            Yxs=[0.5],           # Biomass yield coefficient on solute
            Dt=[4.0E-5],        # Aquious solute diffusion through tank fluid
            Db=[1.0E-5],        # Effective solute diffusion through biofilm
            srcS=[(S,X,Lf,t,z,p) -> 0.0],     # Source of solutes
            
            # --------------- #
            # Tank Parameters #
            # --------------- #
            V=0.1,        # Volume of tank [m³]
            A=1,          # Surface area of biofilm [m²]
            Q=1,        # Flowrate through tank [m³/d]

            # ------------------ #
            # Biofilm Parameters #
            # ------------------ #
            Nz=Nzs[n],          # Number of grid points in biofilm
            Lfo=5.0E-5,     # Biofilm initial thickness [m]
            LL=0.0,         # Boundary layer thickness [m]
        )

        t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver

        # Analytic Result
        Xb=p.rho[1].*Pb[:,end]
        phi = sqrt(mumax.*Xb[1].*Lf[end]^2/
                        (p.Db[1]*KM*p.Yxs[1]))
        Sb_ana = St[end].*cosh.(phi*zm/Lf[end])/cosh(phi)

        error[n]=maximum(abs.(Sb_ana .- Sb'))
    end

    for n in eachindex(Nzs)
        @printf("Nz = %4i, Error = %6.3e \n",Nzs[n], error[n])
    end

    # plot( Nzs, error, xaxis=:log, yaxis=:log)
    # plot!(Nzs, 10.0.*Nzs.^(-2))

    # Compute order (slope of line in log-log space)
    x=log.(Nzs);    xbar = mean(x)
    y=log.(error);  ybar = mean(y)
    order = abs(sum((x.-xbar).*(y.-ybar))/sum((x.-xbar).^2))

    @printf("Order = %6.3f \n",order)

    return order
end

order = test_Diffusion()