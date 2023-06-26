using Parameters

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
