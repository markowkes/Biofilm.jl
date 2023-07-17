using Parameters

@with_kw struct biofilmGrid
    # Grid
    z  :: Vector{Float64}
    zm :: Vector{Float64}
    dz :: Float64
end
