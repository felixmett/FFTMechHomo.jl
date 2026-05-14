@testset "LinearIsotropicElastic" begin
    @testset "constructor" begin
        # correct usage
        mat2 = LinearIsotropicElastic{2}(1.0, 0.3)
        mat3 = LinearIsotropicElastic{3}(1.0, 0.3)
        @test mat2.E ≈ 1.0 && mat3.E ≈ 1.0
        @test mat2.nu ≈ 0.3 && mat3.nu ≈ 0.3
        @test mat2 isa LinearIsotropicElastic{2, Float64}
        @test mat3 isa LinearIsotropicElastic{3, Float64}

        # explicit type
        mat32 = LinearIsotropicElastic{3, Float32}(Float32(1.0), Float32(0.3))
        @test mat32 isa LinearIsotropicElastic{3, Float32}

        # inferred type
        mat32 = LinearIsotropicElastic{3}(Float32(1.0), Float32(0.3))
        @test mat32 isa LinearIsotropicElastic{3, Float32}

        # mixed type
        mat = LinearIsotropicElastic{3}(1, 0.3) # Int + Float
        @test mat isa LinearIsotropicElastic{3, Float64}

        # invalid parameters
        @test_throws ArgumentError LinearIsotropicElastic{3}(1.0, 0.5)   # nu = 0.5
        @test_throws ArgumentError LinearIsotropicElastic{3}(0.0, 0.3)     # E = 0
        @test_throws ArgumentError LinearIsotropicElastic{2}(1.0, -1.0)  # nu = -1
        @test_throws ArgumentError LinearIsotropicElastic{4}(1.0, 0.3)   # invalid dim
    end

    @testset "compute_stress!" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.3)
        mu     = mat.E / (2 * (1 + mat.nu))
        lambda = mat.E * mat.nu / ((1 + mat.nu) * (1 - 2mat.nu))

        # ε₁₁ ≠ 0
        strain = zeros(6, 1, 1, 1)
        stress = zeros(6, 1, 1, 1)
        strain[1, 1, 1, 1] = 1.0
        FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1, 1))

        @test stress[1, 1, 1, 1] ≈ (2mu + lambda) * strain[1, 1, 1, 1] # σ₁₁
        @test stress[2, 1, 1, 1] ≈ lambda * strain[1, 1, 1, 1] # σ₂₂
        @test stress[3, 1, 1, 1] ≈ lambda * strain[1, 1, 1, 1] # σ₃₃
        @test all(stress[4:6, 1, 1, 1] .≈ 0.0) # σ₁₂, σ₁₃, σ₂₃ = 0

        # γ₁₂ ≠ 0
        strain = zeros(6, 1, 1, 1)
        stress = zeros(6, 1, 1, 1)
        strain[4, 1, 1, 1] = 1.0
        FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1, 1))

        @test all(stress[1:3, 1, 1, 1] .≈ 0.0) # σ₁₁, σ₂₂, σ₃₃ = 0
        @test stress[4, 1, 1, 1] ≈ mu * strain[4, 1, 1, 1] # σ₁₂
        @test stress[5, 1, 1, 1] ≈ 0.0 # σ₁₃ = 0
        @test stress[6, 1, 1, 1] ≈ 0.0 # σ₂₃ = 0
    end
end