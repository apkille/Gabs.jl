@testitem "PDMats Benchmarks" begin
    using Gabs
    using LinearAlgebra
    using BenchmarkTools
    using PDMats
    
    function run_benchmarks(nmodes)
        basis = QuadPairBasis(nmodes)
        
        std_state = vacuumstate(basis)
        pd_state = vacuumstate_pd(basis)
        
        testvec = randn(2*nmodes)
        
        std_purity = @benchmark purity($std_state)
        pd_purity = @benchmark purity($pd_state)
        
        std_whiten = @benchmark whiten($std_state, $testvec)
        pd_whiten = @benchmark whiten($pd_state, $testvec)
        
        std_invquad = @benchmark invquad($std_state, $testvec)
        pd_invquad = @benchmark invquad($pd_state, $testvec)
        
        std_create = @benchmark vacuumstate($basis)
        pd_create = @benchmark vacuumstate_pd($basis)
        
        return Dict(
            "purity" => (std_purity, pd_purity),
            "whiten" => (std_whiten, pd_whiten),
            "invquad" => (std_invquad, pd_invquad),
            "creation" => (std_create, pd_create)
        )
    end
    
  
    if get(ENV, "RUN_BENCHMARKS", "false") == "true"
        @testset "Performance comparison" begin
            results_small = run_benchmarks(5)
            results_medium = run_benchmarks(10)
            results_large = run_benchmarks(20)
            
            # for the sake of avoiding test failures
            @test true
            
            # results 
            for (size, results) in zip(["Small (5 modes)", "Medium (10 modes)", "Large (20 modes)"], 
                                      [results_small, results_medium, results_large])
                println("\n===== $size =====")
                for (name, (std, pd)) in results
                    speedup = minimum(std.times) / minimum(pd.times)
                    println("$name: $(round(speedup, digits=2))x speedup")
                end
            end
        end
    else
        @testset "Skip benchmarks" begin
            basis = QuadPairBasis(2)
            std_state = vacuumstate(basis)
            pd_state = vacuumstate_pd(basis)
            
            @test purity(pd_state) ≈ purity(std_state)
            
            testvec = randn(4)
            @test whiten(pd_state, testvec) ≈ whiten(std_state, testvec)
            @test invquad(pd_state, testvec) ≈ invquad(std_state, testvec)
        end
    end
end