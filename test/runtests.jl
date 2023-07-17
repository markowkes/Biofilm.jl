using Biofilm
using Test
using Printf

function testValue(value,expected,tol)
    error = maximum(abs.(value.-expected))
    @printf("Result has error of = %6.4g \n",error)
    return  error <= tol 
end

@testset "Unit Tests" begin

    # Redirect stdout to file
    # open("output_unittests.txt","w") do out
    #     redirect_stdout(out) do 
            println("\n Running test_zeroLL\n ==============")
            include("test_zeroLL.jl")
            println("Checking result") 
            @test computed ≈ analytic atol=1e-12

            println("\n Running test_S_top\n ==============")
            include("test_S_top.jl")
            println("Checking result") 
            @test computed1 ≈ computed2 atol=1e-12

            println("\n Running test_Diffusion\n ==============")
            include("test_Diffusion.jl")
            println("Checking result")
            @test order > 1.5

            println("\n Running test_SteadyState\n ==============")
            include("test_SteadyState.jl")
            println("Checking result") 
            @test computed ≈ analytic atol=1e-2

            println("\n Running test_TimeIntegration\n ==============")
            include("test_TimeIntegration.jl")
            println("Checking result") 
            @test error <= 1e-1
    #     end
    # end

end

@testset "Examples" begin

    # Redirect stdout to file
    # open("output_examples.txt","w") do out
    #    redirect_stdout(out) do 
  
            println("\n Running Case 1\n ==============")
            @testset "Case 1" begin
                include("Case1.jl")
                @test Xt[end] ≈ 256.87   atol=0.1
                @test St[end] ≈ 2.92     atol=0.1
                @test Lf[end] ≈ 0.000309 atol=1e-4
                println(sol.alg)
            end
            
            println("\n Running Case 2\n ==============")
            @testset "Case 2" begin
                include("Case2.jl")
                @test St[1,end] ≈ 2.981401389120606 atol=0.1
                @test St[2,end] ≈ 24.99999970508142 atol=0.1
                @test Lf[end] ≈ 2.998065034446994e-5 atol=1e-4
                println(sol.alg)
            end

            println("\n Running Case 3\n ==============")
            @testset "Case 3" begin
                include("Case3.jl")
                @test Pb[1,end] ≈ 0.072862 atol=1e-3
                @test Pb[2,end] ≈ 0.007137 atol=1e-3
                println(sol.alg)
            end

            println("\n Running Case 4\n ==============")
            @testset "Case 4" begin
                include("Case4.jl")
                @test Xt[1,end] ≈ 6.69736 atol=0.1
                @test Xt[2,end] ≈ 3.80824 atol=0.1
                @test St[1,end] ≈ 10.6042 atol=0.1
                @test St[2,end] ≈ 9.91261 atol=0.1
                println(sol.alg)
            end

            println("\n Running Case 5\n ==============")
            @testset "Case 5" begin
                include("Case5_Phototroph.jl")
                @test Xt[1,end] ≈ 0.230656800294912 atol=0.1
                @test St[1,end] ≈ 9.561670776072377 atol=0.1
                println(sol.alg)
            end

    #     end
    # end
end

@testset "Postprocessing Functions" begin

    # Redirect stdout to file
    # open("output_unittests.txt","w") do out
    #     redirect_stdout(out) do 

            # Run simulation 
            include("Case1.jl")
            
            println("\n Testing biofilm_analyze \n ==============")
            @test_nowarn biofilm_analyze(sol,p,0.5)

            println("\n Testing biofilm_plot \n ==============")
            @test_nowarn biofilm_plot(sol,p)

            println("\n Testing biofilm_analyze \n ==============")
            @test_nowarn biofilm_sol2csv(sol,p; filename="unit_test_123.csv")
            rm("unit_test_123.csv")
            
    #     end
    # end

end
