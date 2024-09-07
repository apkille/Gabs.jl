@testitem "Wigner" begin
    using Gabs
    using StaticArrays

    @testset "symplectic form" begin
        Omega = symplecticform(2)
        Omega_static = symplecticform(SMatrix{4, 4}, 2)
        test_Omega = [0.0 0.0 1.0 0.0;
                      0.0 0.0 0.0 1.0;
                     -1.0 0.0 0.0 0.0;
                      0.0 -1.0 0.0 0.0]
        @test isequal(Omega, test_Omega)
        @test Omega_static isa SMatrix
    end

    @testset "wigner functions" begin
        vac = vacuumstate()
        c = coherentstate(1.0)
        @test isapprox(wignerchar(vac, [0.0, 0.0]), 1.0 - 0.0im)
        @test wigner(vac, [rand(), rand()]) > 0.0
        @test wigner(c, [rand(), rand()]) > 0.0
    end
end