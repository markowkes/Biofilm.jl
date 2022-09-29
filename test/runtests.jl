using Biofilm
using Test

function testValue(value,expected,tol)
    return maximum(abs.(value.-expected)) <= tol 
end

@testset "Biofilm.jl" begin

    # Test Caes 1
    println("\n Running Case 1\n ==============")
    include("Case1.jl")
    println("Checking result")
    @test testValue( X[end],  256.87,0.1) && 
          testValue( S[end],    2.92,0.1) &&
          testValue(Lf[end],0.000309,1e-4)

    println("\n Running Case 2\n ==============")
    include("Case2.jl")
    println("Checking result")
    @test testValue( S[:,end],[2.981401389120606,24.99999970508142],0.1) &&
          testValue(Lf[end],2.998065034446994e-5,1e-4)

    println("\n Running Case 3\n ==============")
    include("Case3.jl")
    println("Checking result")
    @test testValue(Pb[:,1],[0.00632,0.0737],1e-3)

    println("\n Running Case 4\n ==============")
    include("Case4.jl")
    println("Checking result")
    @test testValue( X[:,end],[6.69736,3.80824],0.1) &&
          testValue( S[:,end],[10.6042,9.91261],0.1)

    println("\n Running Case 5\n ==============")
    include("Case5_Phototroph.jl")
    println("Checking result")
    @test testValue( X[end],0.230656800294912,0.1) &&
          testValue( S[end],9.561670776072377,0.1)
end
