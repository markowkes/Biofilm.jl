using Biofilm
using Test

include("Case1.jl")

@testset "Biofilm.jl" begin

    # Test Caes 1
    println("Running Case 1")
    t,X,S,Pb,Sb,Lf=Case1()
    println("Checking result")
    @test test_Case1(t,X,S,Pb,Sb,Lf)

    # Multiple Substrates 
    include("test_multipleSubstrate.jl")
end
