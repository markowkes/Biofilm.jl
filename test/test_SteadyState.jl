# test_SteadyState.jl
# -----------------------------------------------
# Test with large diffusivity with an analytic solution
# -----------------------------------------------

using Biofilm
using Printf

# Loop over # grid points
Nzs = [10,25,50,75,100]

error = similar(Nzs,Float64)
for n in eachindex(Nzs)

    # Constants used for growthrates of particulate(s)
    mumax = 2000;
    KM = 2500;

    # Define a structure to hold all the parameters
    p = param(
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
        Kdet=1900.0,       # Particulates detachment coefficient
        srcX=[(S,X,t,p) -> 0.0],     # Source of particulates
        # Growthrates for each particulate (constants defined above!)
        mu=[(S,X,Lf,t,z,p) -> (mumax * S[1,:]) ./ (KM)],

        # -------------------- #
        # Substrate Parameters #
        # -------------------- #
        SNames=["Oxygen"],   # Substrate names
        Sin=[(t) -> 100],    # Substrate inflow (can be function of time)
        Sto=[25.0],          # Tank substrate concentraion initial condition(s)
        Sbo=[0.0],           # Biofilm substrates concentration initial condition(s)
        Yxs=[0.5],           # Biomass yield coefficient on substrate
        Daq=[400.0],        # Substrate diffusion through boundary layer
        De =[100],        # Substrate diffusion through biofilm     
        srcS=[(S,X,t,p) -> 0.0],     # Source of substrates
        
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
        Lfo=5.0E-6,     # Biofilm initial thickness [m]
        LL=1e-4,         # Boundary layer thickness [m]
    )

    t,zm,Xt,St,Pb,Sb,Lf,sol = BiofilmSolver(p) # Run solver

    # Analytic Result
    begin
        St_ana=[0.0]
        err = 1.0
        while abs(err) > 1e-12
            mu = p.mu[1](St_ana,Xt[:,end],Lf[end],t[end],zm,p)[1]
            Lf = mu / p.Kdet
            Vdet = mu * Lf 
            x = p.Yxs[1]*(p.Sin[1](t[end]) - St_ana[1])
            LHS = p.Q * x 
            RHS = Vdet*p.A*Pb[1,end]*p.rho[1] + mu*x*p.V
            err = LHS - RHS
            St_ana[1] += 0.01*err
        end
    end
    
    # Compare computed and analytic 
    error[n] = abs(St[1,end] - St_ana[1])
    @printf("Computed: %6.3f, Analytic %6.3f, Error = %6.3g \n",St[1,end],St_ana[1],error[n])
    
end


for n in eachindex(Nzs)
    @printf("Nz = %4i, Error = %6.3e \n", Nzs[n], error[n])
end

# using Plots
# plot( Nzs, error, xaxis=:log, yaxis=:log)
# plot!(Nzs, 10.0.*Nzs.^(-2))

# Compute order (slope of line in log-log space)
using Statistics
x=log.(Nzs);    xbar = mean(x)
y=log.(error);  ybar = mean(y)
order = abs(sum((x.-xbar).*(y.-ybar))/sum((x.-xbar).^2))

@printf("Order = %6.3f \n",order)