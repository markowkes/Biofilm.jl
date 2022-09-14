    # Case 1
    
    function Case1()
        # Constants used for growthrates of particulate(s)
        mumax = 20;
        KM = 3;

        # Define a structure to hold all the parameters
        p = param(
            # Growthrates for each particulate (constants defined above!)
            mu=[(S, X, Lf, t, z, p) -> (mumax * S) ./ (KM .+ S)],

            # Source of particulates
            src=[(S, X, p) -> 0.0],

            # Substrate inflow (can be function of time)
            Sin=[(t) -> 100],

            # Time
            tFinal=1,   # Simulation time [days]
            outPeriod=1e-1,  # Time between outputs [days]

            # Simulation
            Title="Single S P Case",
            SNames=["Oxygen"],
            XNames=["Bug"],
            makePlots=false,

            # Tank Geometry
            V=0.1,        # Volume of tank [m³]
            A=1,          # Surface area of biofilm [m²]
            Q=1,          # Flowrate through tank [m³/s]
            Xo=[10.0],     # Tank particulate initial condition(s)
            So=[10.0],    # Tank substrate initial condition(s)
            LL=1.00E-7,    # Boundary layer thickness [m]

            # Biofilm
            Nz=50,            # Number of grid points to represent biofilm
            Pbo=[0.08],     # Biofilm particulates initial condition(s)
            Sbo=[0.0],     # Biofilm substrates initial condition(s)
            Lfo=1.0E-5,    # Biofilm initial thickness [m]

            # Substance Constants
            Yxs=[2.646],     # Biomass yield coeffficient on substrate
            Daq=[4.0E-5],    # Substrate diffusion through boundary layer
            De=[6.9E-5],    # Substrate diffusion through biofilm     
            rho=[2.0E4],     # Particulate densities
            Kdet=20000.0,     # Particulates detachment coefficient

            #Sin_f = 100,
            #Sin_period = 0,         # Substrates concentration(s) into tank

            # Tolerance
            tol=1e-12,
        )

        sol,t,X,S,Pb,Sb,Lf = BiofilmSolver(p) # Run solver
        #outputs(t,X,S,Pb,Sb,Lf,p) # Plot final results

        return t,X,S,Pb,Sb,Lf
    end

    function test_Case1(t,X,S,Pb,Sb,Lf)
        return X[end]≈256.8670111642846 && S[end]≈2.9289825583446394
    end