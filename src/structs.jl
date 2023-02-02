using Parameters

# Define parameter structure with defalut values
@with_kw struct param 
    # Simulation
    tFinal :: Float64         # Simulation time [days]
    outPeriod :: Float64       # Output period [days] 
    tol :: Float64 = 1e-4      # Tolerance
    plotPeriod :: Float64 = outPeriod # Plots period [days]
    makePlots :: Bool = true   # Make plots during simulation 
    discontinuityPeriod :: Float64 = Inf # Period between discontinuites 
    optionalPlot :: String = "growthrate" # What to put in 6th plot
    plotSize :: Tuple{Int64, Int64} = (1600,1000) # Plot size

    # Names
    Title :: String            # Description of case (used on plots)
    SNames :: Vector{String}   # Substrate names (used on plots)
    XNames :: Vector{String}   # Particulate names (used on plots)

    # Grid
    Nz :: Int64               # Number of grid points
    
    # Geometry
    V :: Float64              # Volume of tank
    A :: Float64              # Surface area of biofilm
    LL :: Float64             # Boundary layer thickness
    
    # Flow
    Q :: Float64              # Flowrate through tank
    
    # Initial conditions
    Xto :: Vector{Float64}    # Tank particulate concentrations
    Sto :: Vector{Float64}    # Tank substrate concentrations
    Pbo :: Vector{Float64}    # Biofilm particulate volume fractions
    Sbo :: Vector{Float64}    # Biofilm substrate concentrations
    Lfo :: Float64            # Biofilm thickness
    
    # Substrate parameters
    Yxs :: Array{Float64}  # Biomass yield coefficient on substrate
    Daq :: Vector{Float64}  # Substrate diffusion through boundary layer
    De  :: Vector{Float64}    # Substrate diffusion through biofilm     
    rho :: Vector{Float64}  # Particulate densities
    Kdet :: Float64 # Particulate detachment coefficient
    Sin :: Vector{Function}  # Inflow substrate concentration 
    srcS :: Vector{Function}  # Array of particulate source expressions

    # Particulates
    mu :: Vector{Function}  # Array of growthrate expressions
    srcX :: Vector{Function}  # Array of particulate source expressions
    
    # Computed quantites
    Ptot :: Float64 = sum(Pbo) # Total particulate volume fraction
    Nx :: Int64 = length(Xto)   # Number of particulates
    Ns :: Int64 = length(Sto)   # Number of substrates
    
end

@with_kw struct ranges
    # Ranges - used for spliting dependent variables in ODE solver
    Xt :: UnitRange{Int64}
    St :: UnitRange{Int64} 
    Pb :: UnitRange{Int64}
    Sb :: UnitRange{Int64}
    Lf :: UnitRange{Int64}
end

@with_kw struct biofilmGrid
    # Grid
    z  :: Vector{Float64}
    zm :: Vector{Float64}
    dz :: Float64
end
