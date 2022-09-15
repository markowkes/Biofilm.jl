using Biofilm
using Test

function testValue(value,expected,tol)
    return maximum(abs.(value.-expected)) <= tol 
end

@testset "Biofilm.jl" begin

    # Test Caes 1
    println("\n Running Case 1\n ==============")
    include("../examples/Case1.jl")
    println("Checking result")
    @test testValue( X[end],  256.87,0.1) && 
          testValue( S[end],    2.92,0.1) &&
          testValue(Lf[end],0.000309,1e-4)

    println("\n Running Case 2\n ==============")
    include("../examples/Case2.jl")
    println("Checking result")
    @test testValue( S[:,end],[2.772,25.0],0.1) &&
          testValue(Lf[end],0.00010747,1e-4)

    println("\n Running Case 3\n ==============")
    include("../examples/Case3.jl")
    println("Checking result")
    @test testValue(Pb[:,1],[0.0068849,0.0743802],1e-3)

    println("\n Running Case 4\n ==============")
    include("../examples/Case4.jl")
    println("Checking result")
    @test testValue( X[:,end],[6.69736,3.80824],0.1) &&
          testValue( S[:,end],[10.6042,9.91261],0.1)

    println("\n Running Case 5\n ==============")
    include("../examples/Case5_Phototroph.jl")
    println("Checking result")
    @test testValue( X[end],98.269,0.1) &&
          testValue( S[end],   8.6,0.1)
end
