using Parameters

# Define parameter structure with defalut values
@with_kw struct param 
    # Simulation
    tFinal          # Simulation time [days]
    outPeriod       # Output period [days] 
    tol = 1e-4      # Tolerance
    plotPeriod = outPeriod # Plots period [days]
    makePlots = true   # Make plots during simulation 
    discontinuityPeriod = Inf # Period between discontinuites 

    # Names
    Title            # Description of case (used on plots)
    SNames  # Substrate names (used on plots)
    XNames  # Particulate names (used on plots)

    # Grid
    Nz               # Number of grid points
    
    # Geometry
    V              # Volume of tank
    A              # Surface area of biofilm
    LL             # Boundary layer thickness
    
    # Flow
    Q              # Flowrate through tank
    
    # Initial conditions
    Xo    # Tank particulate concentrations
    So    # Tank substrate concentrations
    Pbo   # Biofilm particulate volume fractions
    Sbo   # Biofilm substrate concentrations
    Lfo   # Biofilm thickness
    
    # Substrate parameters
    Yxs   # Biomass yield coefficient on substrate
    Daq   # Substrate diffusion through boundary layer
    De    # Substrate diffusion through biofilm     
    rho   # Particulate densities
    Kdet  # Particulate detachment coefficient
    Sin   # Inflow substrate concentration 

    
    # Particulates
    mu   # Array of growthrate expressions
    src  # Array of particulate source expressions
    
    # Computed quantites
    Ptot = sum(Pbo) # Total particulate volume fraction
    Nx = length(Xo)   # Number of particulates
    Ns = length(So)   # Number of substrates
end

@with_kw struct ranges
    # Ranges - used for spliting dependent variables in ODE solver
    X 
    S 
    Pb
    Sb
    Lf
end

@with_kw struct biofilmGrid
    # Grid
    z 
    zm
    dz
end
