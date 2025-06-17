@testitem "Phase 3: Integration and Advanced Features" begin
    using Gabs
    using LinearAlgebra
    using StaticArrays

    nmodes1 = 1
    nmodes2 = 2
    qpairbasis1 = QuadPairBasis(nmodes1)
    qpairbasis2 = QuadPairBasis(nmodes2)
    qblockbasis1 = QuadBlockBasis(nmodes1)
    qblockbasis2 = QuadBlockBasis(nmodes2)

    @testset "Gaussian Operations Integration" begin
        @testset "Gaussian Unitary Operations" begin
            for basis in [qpairbasis1, qblockbasis1]
                coh1 = coherentstate(basis, 1.0 + 0.5im)
                coh2 = coherentstate(basis, -1.0 + 0.3im)
                lc = GaussianLinearCombination(basis, [0.6, 0.8], [coh1, coh2])
                
                disp_op = displace(basis, 0.5 - 0.2im)
                lc_displaced = disp_op * lc
                
                @test lc_displaced isa GaussianLinearCombination
                @test length(lc_displaced) == 2
                @test lc_displaced.basis == basis
                @test lc_displaced.ħ == lc.ħ
                @test lc_displaced.coefficients == lc.coefficients  
                
                expected_state1 = disp_op * coh1
                expected_state2 = disp_op * coh2
                @test isapprox(lc_displaced.states[1], expected_state1)
                @test isapprox(lc_displaced.states[2], expected_state2)
                
                squeeze_op = squeeze(basis, 0.3, π/4)
                lc_squeezed = squeeze_op * lc
                
                @test lc_squeezed isa GaussianLinearCombination
                @test lc_squeezed.coefficients == lc.coefficients
                @test isapprox(lc_squeezed.states[1], squeeze_op * coh1)
                @test isapprox(lc_squeezed.states[2], squeeze_op * coh2)
                
                phase_op = phaseshift(basis, π/3)
                lc_phase = phase_op * lc
                
                @test lc_phase isa GaussianLinearCombination
                @test lc_phase.coefficients == lc.coefficients
                @test isapprox(lc_phase.states[1], phase_op * coh1)
                @test isapprox(lc_phase.states[2], phase_op * coh2)
            end
        end
        
        @testset "Gaussian Channel Operations" begin
            for basis in [qpairbasis1, qblockbasis1]
                coh1 = coherentstate(basis, 1.0)
                coh2 = coherentstate(basis, -1.0)
                lc = GaussianLinearCombination(basis, [0.7, 0.3], [coh1, coh2])
                
                att_channel = attenuator(basis, π/6, 2)
                lc_attenuated = att_channel * lc
                
                @test lc_attenuated isa GaussianLinearCombination
                @test length(lc_attenuated) == 2
                @test lc_attenuated.basis == basis
                @test lc_attenuated.ħ == lc.ħ
                @test lc_attenuated.coefficients == lc.coefficients  
                
                @test isapprox(lc_attenuated.states[1], att_channel * coh1)
                @test isapprox(lc_attenuated.states[2], att_channel * coh2)
                
                amp_channel = amplifier(basis, 0.2, 1.5)
                lc_amplified = amp_channel * lc
                
                @test lc_amplified isa GaussianLinearCombination
                @test lc_amplified.coefficients == lc.coefficients
                @test isapprox(lc_amplified.states[1], amp_channel * coh1)
                @test isapprox(lc_amplified.states[2], amp_channel * coh2)
            end
        end
        
        @testset "Operation Compatibility Checks" begin
            lc1 = GaussianLinearCombination(coherentstate(qpairbasis1, 1.0))
            lc2 = GaussianLinearCombination(coherentstate(qblockbasis1, 1.0))
            op_pair = displace(qpairbasis1, 0.5)
            op_block = displace(qblockbasis1, 0.5)
            
            @test_throws ArgumentError op_block * lc1  
            @test_throws ArgumentError op_pair * lc2   
            
            
            lc_h1 = GaussianLinearCombination(coherentstate(qpairbasis1, 1.0, ħ=1))
            op_h2 = displace(qpairbasis1, 0.5, ħ=2)
            
            @test_throws ArgumentError op_h2 * lc_h1  
        end
    end

    @testset "Tensor Products" begin
        @testset "Basic Tensor Products" begin
            coh1 = coherentstate(qpairbasis1, 1.0)
            coh2 = coherentstate(qpairbasis1, -1.0)
            vac1 = vacuumstate(qpairbasis1)
            sq1 = squeezedstate(qpairbasis1, 0.3, π/6)
            
            lc1 = GaussianLinearCombination(qpairbasis1, [0.6, 0.4], [coh1, coh2])
            lc2 = GaussianLinearCombination(qpairbasis1, [0.8, 0.2], [vac1, sq1])
            
            lc_tensor = tensor(lc1, lc2)
            
            @test lc_tensor isa GaussianLinearCombination
            @test length(lc_tensor) == 4  
            @test lc_tensor.basis == qpairbasis2  
            @test lc_tensor.ħ == lc1.ħ
            
            expected_coeffs = [0.6*0.8, 0.6*0.2, 0.4*0.8, 0.4*0.2]
            @test isapprox(sort(lc_tensor.coefficients), sort(expected_coeffs))
            
            expected_states = [
                coh1 ⊗ vac1, coh1 ⊗ sq1,
                coh2 ⊗ vac1, coh2 ⊗ sq1
            ]
            
            for expected_state in expected_states
                found = false
                for actual_state in lc_tensor.states
                    if isapprox(actual_state, expected_state, atol=1e-12)
                        found = true
                        break
                    end
                end
                @test found
            end
        end
        
        @testset "Typed Tensor Products" begin
            coh1 = coherentstate(qpairbasis1, 1.0)
            vac1 = vacuumstate(qpairbasis1)
            lc1 = GaussianLinearCombination(qpairbasis1, [0.7, 0.3], [coh1, vac1])
            lc2 = GaussianLinearCombination(qpairbasis1, [1.0], [coh1])
            
            lc_tensor_typed = tensor(Vector{Float64}, Matrix{Float64}, lc1, lc2)
            @test lc_tensor_typed isa GaussianLinearCombination
            @test length(lc_tensor_typed) == 2  
            @test lc_tensor_typed.coefficients == [0.7, 0.3]
            
            lc_tensor_single = tensor(Matrix{Float64}, lc1, lc2)
            @test lc_tensor_single isa GaussianLinearCombination
            @test length(lc_tensor_single) == 2
        end
        
        @testset "Tensor Product Properties" begin
            vac = vacuumstate(qpairbasis1)
            coh = coherentstate(qpairbasis1, 1.0)
            sq = squeezedstate(qpairbasis1, 0.2, π/4)
            
            lc1 = GaussianLinearCombination(vac)
            lc2 = GaussianLinearCombination(coh)  
            lc3 = GaussianLinearCombination(sq)
            
            left_assoc = tensor(tensor(lc1, lc2), lc3)
            right_assoc = tensor(lc1, tensor(lc2, lc3))
            
            @test left_assoc.basis.nmodes == right_assoc.basis.nmodes == 3
            @test length(left_assoc) == length(right_assoc) == 1
            
            lc_self = tensor(lc1, lc1)
            @test lc_self.basis.nmodes == 2
            @test length(lc_self) == 1
        end
        
        @testset "Tensor Product Compatibility Checks" begin
            lc_pair = GaussianLinearCombination(coherentstate(qpairbasis1, 1.0))
            lc_block = GaussianLinearCombination(coherentstate(qblockbasis1, 1.0))
            lc_h1 = GaussianLinearCombination(coherentstate(qpairbasis1, 1.0, ħ=1))
            lc_h2 = GaussianLinearCombination(coherentstate(qpairbasis1, 1.0, ħ=2))
            
            @test_throws ArgumentError tensor(lc_pair, lc_block)
            
            @test_throws ArgumentError tensor(lc_h1, lc_h2)
        end
        
        @testset "StaticArrays Compatibility" begin
            coh_static = coherentstate(SVector{2}, SMatrix{2,2}, qpairbasis1, 1.0)
            vac_static = vacuumstate(SVector{2}, SMatrix{2,2}, qpairbasis1)
            
            lc1_static = GaussianLinearCombination(qpairbasis1, [0.6, 0.8], [coh_static, vac_static])
            lc2_static = GaussianLinearCombination(qpairbasis1, [1.0], [coh_static])
            
            lc_tensor_static = tensor(lc1_static, lc2_static)
            @test lc_tensor_static isa GaussianLinearCombination
            @test length(lc_tensor_static) == 2
        end
    end

    @testset "Partial Traces" begin
        @testset "Basic Partial Trace" begin
            basis = qpairbasis2
            coh1 = coherentstate(QuadPairBasis(1), 1.0)
            coh2 = coherentstate(QuadPairBasis(1), -1.0)
            vac = vacuumstate(QuadPairBasis(1))
            
            state1 = coh1 ⊗ vac
            state2 = coh2 ⊗ vac  
            state3 = vac ⊗ coh1
            
            lc = GaussianLinearCombination(basis, [0.5, 0.3, 0.2], [state1, state2, state3])
            
            lc_traced = ptrace(lc, 1)
            
            @test lc_traced isa GaussianLinearCombination
            @test lc_traced.basis.nmodes == 1  
            @test lc_traced.ħ == lc.ħ
            
 
            @test length(lc_traced) <= 2  
            
            lc_traced2 = ptrace(lc, 2)
            @test lc_traced2 isa GaussianLinearCombination
            @test lc_traced2.basis.nmodes == 1
            
            @test_throws ArgumentError ptrace(lc, [1, 2])
        end
        
        @testset "Partial Trace with Multiple Indices" begin
            basis3 = QuadPairBasis(3)
            vac = vacuumstate(QuadPairBasis(1))
            coh = coherentstate(QuadPairBasis(1), 1.0)
            
            state = vac ⊗ coh ⊗ vac
            lc = GaussianLinearCombination(state)
            
            lc_traced = ptrace(lc, [1, 3])
            
            @test lc_traced isa GaussianLinearCombination
            @test lc_traced.basis.nmodes == 1
            @test length(lc_traced) == 1
            @test isapprox(lc_traced.states[1], coh, atol=1e-12)
        end
        
        @testset "Typed Partial Trace" begin
            basis = qpairbasis2
            state1 = coherentstate(QuadPairBasis(1), 1.0) ⊗ vacuumstate(QuadPairBasis(1))
            lc = GaussianLinearCombination(basis, [1.0], [state1])
            
            lc_traced_typed = ptrace(Vector{Float64}, Matrix{Float64}, lc, 1)
            @test lc_traced_typed isa GaussianLinearCombination
            @test lc_traced_typed.basis.nmodes == 1
            
            lc_traced_single = ptrace(Matrix{Float64}, lc, 1)
            @test lc_traced_single isa GaussianLinearCombination
        end
        
        @testset "Partial Trace Error Handling" begin
            basis = qpairbasis2
            coh = coherentstate(QuadPairBasis(1), 1.0)
            vac = vacuumstate(QuadPairBasis(1))
            state = coh ⊗ vac
            lc = GaussianLinearCombination(basis, [1.0], [state])
            
            @test_throws ArgumentError ptrace(lc, 3)  # ..Index too large
            @test_throws ArgumentError ptrace(lc, 0)  # Index too small
            @test_throws ArgumentError ptrace(lc, [1, 2])  # Tracing all modes
            @test_throws ArgumentError ptrace(lc, [1, 1])  # Duplicate indices
        end
        
        @testset "Partial Trace Simplification" begin
            basis = qpairbasis2
            coh = coherentstate(QuadPairBasis(1), 1.0)
            vac = vacuumstate(QuadPairBasis(1))
            
            state1 = coh ⊗ vac
            state2 = coh ⊗ squeezedstate(QuadPairBasis(1), 0.1, 0.0) 
            
            lc = GaussianLinearCombination(basis, [0.6, 0.4], [state1, state2])
            
            # note - Trace out the second mode - both should trace to |coh⟩
            lc_traced = ptrace(lc, 2)
            
            @test lc_traced isa GaussianLinearCombination
     
            @test lc_traced.basis.nmodes == 1
        end
    end

    @testset "Wigner Functions with Quantum Interference" begin
        @testset "Cross-Wigner Functions" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            x = [0.5, 0.3]  
            
            cross_w = cross_wigner(coh1, coh2, x)
            @test cross_w isa ComplexF64
            @test isfinite(cross_w)
            

            cross_w_reverse = cross_wigner(coh2, coh1, x)
            @test isapprox(cross_w, conj(cross_w_reverse), rtol=1e-8)
            
           
            cross_w_same = cross_wigner(coh1, coh1, x)
            regular_w = wigner(coh1, x)
            
            identity_holds = isapprox(real(cross_w_same), regular_w, rtol=1e-6)
            @test identity_holds || begin
                println("WARNING: Cross-Wigner identity test failed!")
                println("cross_wigner(s,s,x) = ", cross_w_same)
                println("wigner(s,x) = ", regular_w)
                println("Ratio = ", real(cross_w_same)/regular_w)
                println("This indicates a potential bug in the cross_wigner implementation")
                false 
            end
            
            @test abs(imag(cross_w_same)) < 1e-10
            
            vac = vacuumstate(basis)
            sq = squeezedstate(basis, 0.3, π/4)
            
            cross_w_vac_sq = cross_wigner(vac, sq, x)
            @test cross_w_vac_sq isa ComplexF64
            @test isfinite(cross_w_vac_sq)
        end
        
        @testset "Wigner Function with Interference" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            cat_state = GaussianLinearCombination(basis, [0.7, 0.3], [coh1, coh2])
            
            x = [0.0, 0.0] 
            
            w_interference = wigner(cat_state, x)
            @test w_interference isa Real
            @test isfinite(w_interference)
            
            w1 = wigner(coh1, x)
            w2 = wigner(coh2, x)
            w_cross = cross_wigner(coh1, coh2, x)
            
            expected_w = abs2(0.7) * w1 + abs2(0.3) * w2 + 2 * real(conj(0.7) * 0.3 * w_cross)
            
            @test isapprox(w_interference, expected_w, atol=1e-12)
            

            symmetric_cat = GaussianLinearCombination(basis, [0.5, 0.5], [coh1, coh2])
            x_interference = [0.0, 0.5] 
            w_symmetric = wigner(symmetric_cat, x_interference)
            
            @test isfinite(w_symmetric)
        end
        
        @testset "Cross-Wigner Characteristic Functions" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0 + 0.5im)
            coh2 = coherentstate(basis, -0.5 + 1.0im)
            xi = [0.2, 0.3]  
            
            cross_char = cross_wignerchar(coh1, coh2, xi)
            @test cross_char isa ComplexF64
            @test isfinite(cross_char)
            
       
            cross_char_reverse = cross_wignerchar(coh2, coh1, -xi) 
            hermitian_holds = isapprox(cross_char, conj(cross_char_reverse), rtol=1e-8)
            @test hermitian_holds || begin
                println("WARNING: Cross-Wigner characteristic function Hermitian property failed!")
                println("χ₁₂(ξ) = ", cross_char)
                println("χ₂₁*(-ξ) = ", conj(cross_char_reverse)) 
                println("Difference = ", cross_char - conj(cross_char_reverse))
                println("This indicates a potential bug in cross_wignerchar implementation")
                false 
            end
            
            cross_char_same = cross_wignerchar(coh1, coh1, xi)
            regular_char = wignerchar(coh1, xi)
            identity_holds = isapprox(cross_char_same, regular_char, rtol=1e-8)
            @test identity_holds || begin
                println("WARNING: Cross-Wigner characteristic function identity failed!")
                println("cross_wignerchar(s,s,ξ) = ", cross_char_same)
                println("wignerchar(s,ξ) = ", regular_char)
                println("This indicates a potential bug in cross_wignerchar implementation")
                false 
            end
        end
        
        @testset "Wigner Characteristic Function with Interference" begin
            basis = qpairbasis1
            
            vac = vacuumstate(basis)
            coh = coherentstate(basis, 1.0)
            superpos = GaussianLinearCombination(basis, [0.6, 0.8], [vac, coh])
            
            xi = [0.1, 0.2]
            
            char_interference = wignerchar(superpos, xi)
            @test char_interference isa ComplexF64
            @test isfinite(char_interference)
            
            char1 = wignerchar(vac, xi)
            char2 = wignerchar(coh, xi)
            char_cross = cross_wignerchar(vac, coh, xi)
            
            expected_char = abs2(0.6) * char1 + abs2(0.8) * char2 + 
                           2 * real(conj(0.6) * 0.8 * char_cross)
                           
            @test isapprox(real(char_interference), real(expected_char), rtol=1e-10)
        end
        
        @testset "Wigner Function Error Handling" begin
            basis = qpairbasis1
            coh = coherentstate(basis, 1.0)
            lc = GaussianLinearCombination(coh)
            
            @test_throws ArgumentError wigner(lc, [1.0])  
            @test_throws ArgumentError wigner(lc, [1.0, 2.0, 3.0])  
            @test_throws ArgumentError wignerchar(lc, [1.0])
            
            cross_w = cross_wigner(coh, coh, [0.5, 0.3])
            @test cross_w isa ComplexF64
            @test isfinite(cross_w)
            
            @test_throws ArgumentError wignerchar(lc, [1.0, 2.0, 3.0])  
        end
        
        @testset "Multi-mode Wigner Functions" begin
            basis = qpairbasis2
            coh1_2mode = coherentstate(QuadPairBasis(1), 1.0) ⊗ vacuumstate(QuadPairBasis(1))
            coh2_2mode = coherentstate(QuadPairBasis(1), -1.0) ⊗ vacuumstate(QuadPairBasis(1))
            
            lc_2mode = GaussianLinearCombination(basis, [0.7, 0.3], [coh1_2mode, coh2_2mode])
            
            x_2mode = [0.1, 0.2, 0.3, 0.4] 
            
            w_2mode = wigner(lc_2mode, x_2mode)
            @test w_2mode isa Real
            @test isfinite(w_2mode)
            
            xi_2mode = [0.05, 0.1, 0.15, 0.2]
            char_2mode = wignerchar(lc_2mode, xi_2mode)
            @test char_2mode isa ComplexF64
            @test isfinite(char_2mode)
        end
    end

    @testset "Advanced State Metrics" begin
        @testset "Purity of Linear Combinations" begin
            basis = qpairbasis1
            
            coh = coherentstate(basis, 1.0)
            lc_single = GaussianLinearCombination(coh)
            @test purity(lc_single) == 1.0
            
            vac = vacuumstate(basis)
            lc_super = GaussianLinearCombination(basis, [0.6, 0.8], [coh, vac])
            @test purity(lc_super) == 1.0  
            
            cat_state = catstate_even(basis, 1.0)
            @test purity(cat_state) == 1.0
            
            gkp_state = gkpstate(basis, lattice=:square, delta=0.2, nmax=2)
            @test purity(gkp_state) == 1.0
            
            lc_complex = GaussianLinearCombination(basis, [0.6 + 0.8im, 0.2 - 0.1im], [coh, vac])
            @test purity(lc_complex) == 1.0
        end
        
        @testset "Von Neumann Entropy of Linear Combinations" begin
            basis = qpairbasis1
            
            sq = squeezedstate(basis, 0.3, π/6)
            lc_single = GaussianLinearCombination(sq)
            @test entropy_vn(lc_single) == 0.0
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            lc_super = GaussianLinearCombination(basis, [0.7, 0.3], [coh1, coh2])
            @test entropy_vn(lc_super) == 0.0  
            
            cat_odd = catstate_odd(basis, 1.5)
            @test entropy_vn(cat_odd) == 0.0
            
            states = [coherentstate(basis, Float64(i)) for i in 1:5]
            coeffs = [1.0/√5 for _ in 1:5] 
            lc_many = GaussianLinearCombination(basis, coeffs, states)
            @test entropy_vn(lc_many) == 0.0
            
            lc_complex = GaussianLinearCombination(basis, [0.6im, 0.8], [coh1, coh2])
            @test entropy_vn(lc_complex) == 0.0
        end
    end

    @testset "Measurement Theory" begin
        @testset "Basic Measurement Probability" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            lc = GaussianLinearCombination(basis, [0.6, 0.8], [coh1, coh2])
            
            measurement_state = coherentstate(basis, 0.5)
            
            prob = measurement_probability(lc, measurement_state, 1)
            
            @test prob isa Real
            @test 0 <= prob <= 1 
            @test isfinite(prob)
            
 
            norm_squared = abs2(0.6) + abs2(0.8) 
            normalized_lc = GaussianLinearCombination(basis, [0.6, 0.8] ./ sqrt(norm_squared), [coh1, coh2])
            
            prob_self = measurement_probability(normalized_lc, coh1, 1)
            @test prob_self >= 0.0  
        end
        
        @testset "Multi-mode Measurement" begin
            basis = qpairbasis2
            
            coh = coherentstate(QuadPairBasis(1), 1.0)
            vac = vacuumstate(QuadPairBasis(1))
            
            state1 = coh ⊗ vac
            state2 = vac ⊗ coh
            
            lc = GaussianLinearCombination(basis, [0.7, 0.3], [state1, state2])
            
            measurement_1mode = coherentstate(QuadPairBasis(1), 0.5)
            prob1 = measurement_probability(lc, measurement_1mode, 1)
            
            @test prob1 isa Real
            @test 0 <= prob1 <= 1
            @test isfinite(prob1)
            
            prob2 = measurement_probability(lc, measurement_1mode, 2)
            
            @test prob2 isa Real
            @test 0 <= prob2 <= 1
            @test isfinite(prob2)
            

        end
        
        @testset "Measurement with Partial Trace" begin
            basis = qpairbasis2
            
            coh = coherentstate(QuadPairBasis(1), 1.0)
            vac = vacuumstate(QuadPairBasis(1))
            sq = squeezedstate(QuadPairBasis(1), 0.2, π/4)
            
            state1 = coh ⊗ vac
            state2 = vac ⊗ sq
            
            lc = GaussianLinearCombination(basis, [0.6, 0.8], [state1, state2])
            
            measurement = vacuumstate(QuadPairBasis(1))
            prob = measurement_probability(lc, measurement, [1])
            
            @test prob isa Real
            @test 0 <= prob <= 1
            @test isfinite(prob)
        end
        
        @testset "Measurement Error Handling" begin
            basis = qpairbasis2
            coh = coherentstate(QuadPairBasis(1), 1.0)
            vac = vacuumstate(QuadPairBasis(1))
            state = coh ⊗ vac
            lc = GaussianLinearCombination(basis, [1.0], [state])
            
            measurement_wrong = coherentstate(qpairbasis2, 1.0)  
            @test_throws ArgumentError measurement_probability(lc, measurement_wrong, 1)
            
            measurement_h1 = coherentstate(QuadPairBasis(1), 1.0, ħ=1)
            @test_throws ArgumentError measurement_probability(lc, measurement_h1, 1)
            
            measurement_1mode = coherentstate(QuadPairBasis(1), 1.0)
            @test_throws ArgumentError measurement_probability(lc, measurement_1mode, 3)  
        end
        
        @testset "Born Rule Verification" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 2.0)  
            coh2 = coherentstate(basis, -2.0)
            
            lc = GaussianLinearCombination(basis, [1.0, 1.0], [coh1, coh2])
            
            prob1 = measurement_probability(lc, coh1, 1)
            
            prob2 = measurement_probability(lc, coh2, 1)
            
            @test prob1 >= 0
            @test prob2 >= 0
            @test isfinite(prob1)
            @test isfinite(prob2)
            
   
        end
    end

    @testset "Coherence Measure" begin
        @testset "Basic Coherence Properties" begin
            basis = qpairbasis1
            
            coh = coherentstate(basis, 1.0)
            lc_single = GaussianLinearCombination(coh)
            coherence_single = coherence_measure(lc_single)
            
            @test coherence_single == 1.0
            
            vac = vacuumstate(basis)
            coh2 = coherentstate(basis, 0.1) 
            lc_overlap = GaussianLinearCombination(basis, [0.7, 0.3], [vac, coh2])
            coherence_overlap = coherence_measure(lc_overlap)
            
            @test 0 <= coherence_overlap <= 1
            @test isfinite(coherence_overlap)
        end
        
        @testset "Coherence with Different Overlaps" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 3.0)
            coh2 = coherentstate(basis, -3.0)
            lc_separated = GaussianLinearCombination(basis, [0.7, 0.3], [coh1, coh2])
            coherence_separated = coherence_measure(lc_separated)
            
            coh3 = coherentstate(basis, 0.1)
            coh4 = coherentstate(basis, 0.2)
            lc_close = GaussianLinearCombination(basis, [0.7, 0.3], [coh3, coh4])
            coherence_close = coherence_measure(lc_close)
            
            @test 0 <= coherence_separated <= 1
            @test 0 <= coherence_close <= 1
            
            @test coherence_separated >= coherence_close - 0.1  
        end
        
        @testset "Coherence with Identical States" begin
            basis = qpairbasis1
            
            coh = coherentstate(basis, 1.0)
            lc_identical = GaussianLinearCombination(basis, [0.5, 0.3, 0.2], [coh, coh, coh])
            coherence_identical = coherence_measure(lc_identical)
            
            @test 0 <= coherence_identical <= 1
            @test isfinite(coherence_identical)
            
            @test coherence_identical < 1.0
        end
        
        @testset "Coherence with Many States" begin
            basis = qpairbasis1
            
            states = [coherentstate(basis, Float64(i)) for i in 1:5]
            coeffs = [1.0/5 for _ in 1:5] 
            lc_many = GaussianLinearCombination(basis, coeffs, states)
            coherence_many = coherence_measure(lc_many)
            
            @test 0 <= coherence_many <= 1
            @test isfinite(coherence_many)
            
            @test coherence_many > 0.5 
        end
        
        @testset "Coherence Edge Cases" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            lc_small = GaussianLinearCombination(basis, [1e-10, 1.0], [coh1, coh2])
            coherence_small = coherence_measure(lc_small)
            
            @test 0 <= coherence_small <= 1
            @test isfinite(coherence_small)
            
            many_states = [coherentstate(basis, 0.5 * i) for i in 1:10]
            many_coeffs = [1.0/10 for _ in 1:10]
            lc_many_small = GaussianLinearCombination(basis, many_coeffs, many_states)
            coherence_many_small = coherence_measure(lc_many_small)
            
            @test 0 <= coherence_many_small <= 1
            @test isfinite(coherence_many_small)
        end
    end

    @testset "Integration with Existing Cat and GKP States" begin
        @testset "Cat State Integration" begin
            basis = qpairbasis1
            
            cat_even = catstate_even(basis, 1.0)
            cat_odd = catstate_odd(basis, 1.5)
            
            disp_op = displace(basis, 0.5)
            cat_displaced = disp_op * cat_even
            
            @test cat_displaced isa GaussianLinearCombination
            @test length(cat_displaced) == 2
            
            cat_tensor = tensor(cat_even, cat_odd)
            @test cat_tensor isa GaussianLinearCombination
            @test length(cat_tensor) == 4  
            @test cat_tensor.basis.nmodes == 2
            
            x = [0.0, 0.0]
            w_cat = wigner(cat_even, x)
            @test w_cat isa Real
            @test isfinite(w_cat)
            

        end
        
        @testset "GKP State Integration" begin
            basis = qpairbasis1
            
            gkp_square = gkpstate(basis, lattice=:square, delta=0.2, nmax=2)
            gkp_hex = gkpstate(basis, lattice=:hexagonal, delta=0.2, nmax=1)
            
            
            squeeze_op = squeeze(basis, 0.1, π/6)
            gkp_squeezed = squeeze_op * gkp_square
            
            @test gkp_squeezed isa GaussianLinearCombination
            @test length(gkp_squeezed) == length(gkp_square)
            

            vac = vacuumstate(QuadPairBasis(1))
            gkp_2mode = tensor(gkp_square, GaussianLinearCombination(vac))
            
            @test gkp_2mode.basis.nmodes == 2
            
            gkp_traced = ptrace(gkp_2mode, 2)
            @test gkp_traced.basis.nmodes == 1
            
            measurement = vacuumstate(basis)
            prob = measurement_probability(gkp_square, measurement, 1)
            
            @test 0 <= prob <= 1
            @test isfinite(prob)
        end
        
        @testset "Mixed Cat and GKP Operations" begin
            basis = qpairbasis1
            
            cat = catstate_even(basis, 1.0)
            gkp = gkpstate(basis, lattice=:square, delta=0.3, nmax=1)
            
            mixed = 0.6 * cat + 0.4 * gkp
            @test mixed isa GaussianLinearCombination
            @test mixed.basis == basis
            
            cat_gkp_tensor = tensor(cat, gkp)
            @test cat_gkp_tensor isa GaussianLinearCombination
            @test cat_gkp_tensor.basis.nmodes == 2
            @test length(cat_gkp_tensor) == length(cat) * length(gkp)
            
            x = [0.1, 0.2, 0.3, 0.4]  
            w_mixed = wigner(cat_gkp_tensor, x)
            @test w_mixed isa Real
            @test isfinite(w_mixed)
            
            coherence_mixed = coherence_measure(mixed)
            @test 0 <= coherence_mixed <= 1
        end
    end

    @testset "Performance and Numerical Stability" begin
        @testset "Large Linear Combinations" begin
            basis = qpairbasis1
            
            n_states = 50
            states = [coherentstate(basis, 0.1 * i) for i in 1:n_states]
            coeffs = [1.0/√n_states for _ in 1:n_states]  # Normalized
            
            lc_large = GaussianLinearCombination(basis, coeffs, states)
            
            @test length(lc_large) == n_states
            
            x = [0.5, 0.3]
            w_large = wigner(lc_large, x)
            @test isfinite(w_large)
            
            coherence_large = coherence_measure(lc_large)
            @test 0 <= coherence_large <= 1
            @test isfinite(coherence_large)
            
            @test purity(lc_large) == 1.0
            @test entropy_vn(lc_large) == 0.0
        end
        
        @testset "Extreme Parameter Values" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            lc_tiny = GaussianLinearCombination(basis, [1e-15, 1.0], [coh1, coh2])
            
            @test isfinite(purity(lc_tiny))
            @test isfinite(entropy_vn(lc_tiny))
            
            lc_huge = GaussianLinearCombination(basis, [1e10, 1e10], [coh1, coh2])
            
            x = [0.0, 0.0]
            w_huge = wigner(lc_huge, x)
            @test isfinite(w_huge)
            
            lc_complex = GaussianLinearCombination(basis, [1.0 + 100.0im, 1.0 - 100.0im], [coh1, coh2])
            
            @test isfinite(purity(lc_complex))
            @test entropy_vn(lc_complex) == 0.0
        end
        
        @testset "Numerical Precision" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, 1.0 + 1e-12)  
            
            lc_close = GaussianLinearCombination(basis, [0.7, 0.3], [coh1, coh2])
            
            x = [0.5, 0.5]
            cross_w = cross_wigner(coh1, coh2, x)
            @test isfinite(cross_w)
            @test abs(cross_w) > 0  
            
            w_close = wigner(lc_close, x)
            @test isfinite(w_close)
            
            coherence_close = coherence_measure(lc_close)
            @test 0 <= coherence_close <= 1
            @test isfinite(coherence_close)
        end
    end

    @testset "Mathematical Properties Verification" begin
        @testset "Tensor Product Linearity" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 1.0)
            coh2 = coherentstate(basis, -1.0)
            vac = vacuumstate(basis)
            
            lc1 = GaussianLinearCombination(basis, [0.6, 0.4], [coh1, coh2])
            lc2 = GaussianLinearCombination(basis, [0.8, 0.2], [coh2, vac])
            lc3 = GaussianLinearCombination(vac)
            
            a, b = 0.7, 0.3
            
            left_side = tensor(a * lc1 + b * lc2, lc3)
            right_side = a * tensor(lc1, lc3) + b * tensor(lc2, lc3)
            
            Gabs.simplify!(left_side)
            Gabs.simplify!(right_side)
            
            @test length(left_side) == length(right_side)
            @test isapprox(sum(abs2, left_side.coefficients), sum(abs2, right_side.coefficients), rtol=1e-10)
        end
        
        @testset "Measurement Probability Conservation" begin
            basis = qpairbasis1
            
            coh1 = coherentstate(basis, 3.0)  
            coh2 = coherentstate(basis, -3.0)
            
            lc = GaussianLinearCombination(basis, [1.0, 1.0], [coh1, coh2])
            Gabs.normalize!(lc)
            
            prob1 = measurement_probability(lc, coh1, 1)
            prob2 = measurement_probability(lc, coh2, 1)
            
            @test prob1 >= 0
            @test prob2 >= 0
            @test isfinite(prob1)
            @test isfinite(prob2)
            

            @test prob1 > 0.1  
            @test prob2 > 0.1  
            
            coh_ortho = coherentstate(basis, 3.0im) 
            prob_ortho = measurement_probability(lc, coh_ortho, 1)
            @test prob_ortho >= 0
            @test prob_ortho < 0.5 
        end
        
        @testset "Wigner Function Integration Properties" begin
            basis = qpairbasis1
            
            coh = coherentstate(basis, 0.5)
            vac = vacuumstate(basis)
            lc = GaussianLinearCombination(basis, [0.8, 0.6], [coh, vac])
            
            x_vals = [[0.0, 0.0], [0.5, 0.0], [0.0, 0.5], [1.0, 1.0]]
            
            for x in x_vals
                w_val = wigner(lc, x)
                @test isfinite(w_val)
                
 
                @test abs(w_val) < 1e10 
            end
        end
        
        @testset "Cross-Wigner Symmetry Properties" begin
            basis = qpairbasis1
            
            states = [
                vacuumstate(basis),
                coherentstate(basis, 1.0),
                coherentstate(basis, -1.0),
                squeezedstate(basis, 0.2, π/4)
            ]
            
            x = [0.5, 0.3]
            
            for i in 1:length(states)
                for j in i:length(states)
                    s1, s2 = states[i], states[j]
                    
                    w12 = cross_wigner(s1, s2, x)
                    w21 = cross_wigner(s2, s1, x)
                    
                    @test isapprox(w12, conj(w21), rtol=1e-8)
                    
                    if i == j
                        w_regular = wigner(s1, x)
                        @test isapprox(real(w12), w_regular, rtol=1e-10)
                        @test abs(imag(w12)) < 1e-12  
                    end
                end
            end
        end
        
        @testset "Partial Trace Dimension Consistency" begin
            basis1 = QuadPairBasis(1)
            basis3 = QuadPairBasis(3)
            
            coh = coherentstate(basis1, 1.0)
            vac = vacuumstate(basis1)
            sq = squeezedstate(basis1, 0.3, π/6)
            
            state_3mode = tensor(tensor(coh, vac), sq)
            lc_3mode = GaussianLinearCombination(state_3mode)
            
          
            lc_trace_1 = ptrace(lc_3mode, [1])     
            lc_trace_2 = ptrace(lc_3mode, [2])      
            lc_trace_12 = ptrace(lc_3mode, [1,2])   
            
            @test lc_trace_1.basis.nmodes == 2
            @test lc_trace_2.basis.nmodes == 2
            @test lc_trace_12.basis.nmodes == 1
            
            @test isapprox(lc_trace_12.states[1], sq, rtol=1e-10)
        end
    end
end