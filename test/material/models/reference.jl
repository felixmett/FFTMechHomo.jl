@testset "ReferenceMaterial" begin
    @testset "Constructor" begin
        for dim in (2,3)
            @testset "$(dim)D" begin
                mat = FFTMechHomo.ReferenceMaterial{dim}(1.0)
                @test mat.α₀ ≈ 1.0
                @test mat isa FFTMechHomo.ReferenceMaterial{dim, Float64}

                # explicit type
                mat = FFTMechHomo.ReferenceMaterial{dim, Float32}(Float32(1.0))
                @test mat isa FFTMechHomo.ReferenceMaterial{dim, Float32}
                @test mat.α₀ isa Float32

                # inferred type
                mat = FFTMechHomo.ReferenceMaterial{dim}(Float32(1.0))
                @test mat isa FFTMechHomo.ReferenceMaterial{dim, Float32}

                # from Int
                mat = FFTMechHomo.ReferenceMaterial{dim}(1) # Int
                @test mat isa FFTMechHomo.ReferenceMaterial{dim, Float64}
            end
        end
        
        # invalid parameters
        @test_throws ArgumentError FFTMechHomo.ReferenceMaterial{3}(0.0) # α₀ = 0
        @test_throws ArgumentError FFTMechHomo.ReferenceMaterial{2}(-0.1) # α₀ < 0
        @test_throws ArgumentError FFTMechHomo.ReferenceMaterial{4}(1.0) # invalid dim
    end

    @testset "compute_stress!" begin
        @testset "2D" begin
            mat = FFTMechHomo.ReferenceMaterial{2}(1.0)
            μ = 0.5 * mat.α₀

            # ε₁₁ ≠ 0
            strain = zeros(3, 1, 1)
            stress = zeros(3, 1, 1)
            strain[1, 1, 1] = 1.0
            FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1))

            @test stress[1, 1, 1] ≈ 2μ * strain[1, 1, 1] # σ₁₁
            @test stress[2, 1, 1] ≈ 0.0 # σ₂₂ = 0
            @test stress[3, 1, 1] ≈ 0.0 # σ₁₂ = 0

            # γ₁₂ ≠ 0
            strain = zeros(3, 1, 1)
            stress = zeros(3, 1, 1)
            strain[3, 1, 1] = 1.0
            FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1))

            @test all(stress[1:2, 1, 1] .≈ 0.0) # σ₁₁, σ₂₂ = 0
            @test stress[3, 1, 1] ≈ μ * strain[3, 1, 1] # σ₁₂
        end

        @testset "3D" begin
            mat = FFTMechHomo.ReferenceMaterial{3}(1.0)
            μ = 0.5 * mat.α₀

            # ε₁₁ ≠ 0
            strain = zeros(6, 1, 1, 1)
            stress = zeros(6, 1, 1, 1)
            strain[1, 1, 1, 1] = 1.0
            FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1, 1))

            @test stress[1, 1, 1, 1] ≈ 2μ * strain[1, 1, 1, 1] # σ₁₁
            @test stress[2, 1, 1, 1] ≈ 0.0 # σ₂₂ = 0
            @test stress[3, 1, 1, 1] ≈ 0.0 # σ₃₃ = 0
            @test all(stress[4:6, 1, 1, 1] .≈ 0.0) # σ₂₃, σ₁₃, σ₁₂ = 0

            # γ₂₃ ≠ 0
            strain = zeros(6, 1, 1, 1)
            stress = zeros(6, 1, 1, 1)
            strain[4, 1, 1, 1] = 1.0
            FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1, 1))

            @test all(stress[1:3, 1, 1, 1] .≈ 0.0) # σ₁₁, σ₂₂, σ₃₃ = 0
            @test stress[4, 1, 1, 1] ≈ μ * strain[4, 1, 1, 1] # σ₂₃
            @test stress[5, 1, 1, 1] ≈ 0.0 # σ₁₃ = 0
            @test stress[6, 1, 1, 1] ≈ 0.0 # σ₁₂ = 0
        end
    end
end