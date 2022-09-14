using Parameters

# Define parameter structure with defalut values
@with_kw struct param 
    # Simulation
    tFinal::Float64          # Simulation time [days]
    outPeriod::Float64       # Output period [days] 
    tol::Float64 = 1e-4      # Tolerance
    makePlots::Bool = true   # Make plots during simulation 

    # Names
    Title::String            # Description of case (used on plots)
    SNames::Array{String,1}  # Substrate names (used on plots)
    XNames::Array{String,1}  # Particulate names (used on plots)

    # Grid
    Nz::Int64               # Number of grid points
    
    # Geometry
    V::Float64              # Volume of tank
    A::Float64              # Surface area of biofilm
    LL::Float64             # Boundary layer thickness
    
    # Flow
    Q::Float64              # Flowrate through tank
    
    # Initial conditions
    Xo ::Array{Float64,1}   # Tank particulate concentrations
    So ::Array{Float64,1}   # Tank substrate concentrations
    Pbo::Array{Float64,1}   # Biofilm particulate volume fractions
    Sbo::Array{Float64,1}   # Biofilm substrate concentrations
    Lfo::Float64            # Biofilm thickness
    
    # Substrate parameters
    Yxs ::Array{Float64,2}  # Biomass yield coeffficient on substrate
    Daq ::Array{Float64,1}  # Substrate diffusion through boundary layer
    De  ::Array{Float64,1}  # Substrate diffusion through biofilm     
    rho ::Array{Float64,1}  # Particulate densities
    Kdet::Float64           # Particulate detachment coefficient
    Sin::Array{Function,1}  # Inflow substrate concentration 

    
    # Particulates
    mu ::Array{Function,1}  # Array of growthrate expressions
    src::Array{Function,1}  # Array of particulate source expressions
    
    # Computed quantites
    Ptot::Float64 = sum(Pbo) # Total particulate volume fraction
    Nx::Int64 = length(Xo)   # Number of particulates
    Ns::Int64 = length(So)   # Number of substrates
end

@with_kw struct ranges
    # Ranges - used for spliting dependent variables in ODE solver
    X ::UnitRange{Int64}
    S ::UnitRange{Int64}
    Pb::UnitRange{Int64}
    Sb::UnitRange{Int64}
    Lf::UnitRange{Int64}
end

@with_kw struct biofilmGrid
    # Grid
    z  ::StepRangeLen
    zm ::StepRangeLen
    dz ::Float64
end
