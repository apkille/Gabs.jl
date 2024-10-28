@testitem "Measurements" begin
    using Gabs
    using StaticArrays

    @testset "generaldyne" begin    
        s = vacuumstate()
        state = s ⊗ s ⊗ s ⊗ s
        gd1 = Generaldyne(state, s ⊗ s, [2, 4])
        out1 = outcome(gd1)
        @test isequal(out1, s ⊗ s ⊗ s)

        coh = coherentstate(1.0+im)
        cohs = coh ⊗ s ⊗ coh ⊗ s
        epr = eprstate(1.0, 3.0)
        gd2 = Generaldyne(cohs, epr, [1, 4])
        out2 = outcome(gd2)
        @test isequal(out2, vac ⊗ coh)

        sstatic = vacuumstate(SVector{2}, SMatrix{2,2})
        statestatic = sstatic ⊗ sstatic ⊗ sstatic ⊗ sstatic
        gdstatic = Generaldyne(statestatic, sstatic, [2])
        outstatic = outcome(gdstatic)
        @test isequal(outstatic, sstatic ⊗ sstatic ⊗ sstatic)
    end
end