@testset "LinearIsotropicElastic" begin
    @testset "Constructor" begin
        @testset "2D" begin
            mat = LinearIsotropicElastic{2}(1.0, 0.3)
            @test mat.E ≈ 1.0 && mat.nu ≈ 0.3
            @test mat isa LinearIsotropicElastic{2, Float64}

            # explicit type
            mat = LinearIsotropicElastic{2, Float32}(Float32(1.0), Float32(0.3))
            @test mat isa LinearIsotropicElastic{2, Float32}
            @test mat.E isa Float32
            @test mat.nu isa Float32

            # inferred type
            mat = LinearIsotropicElastic{2}(Float32(1.0), Float32(0.3))
            @test mat isa LinearIsotropicElastic{2, Float32}

            # mixed type
            mat = LinearIsotropicElastic{2}(1, 0.3) # Int + Float
            @test mat isa LinearIsotropicElastic{2, Float64}
        end

        @testset "3D" begin
            mat = LinearIsotropicElastic{3}(1.0, 0.3)
            @test mat.E ≈ 1.0 && mat.nu ≈ 0.3
            @test mat isa LinearIsotropicElastic{3, Float64}

            # explicit type
            mat = LinearIsotropicElastic{3, Float32}(Float32(1.0), Float32(0.3))
            @test mat isa LinearIsotropicElastic{3, Float32}
            @test mat.E isa Float32
            @test mat.nu isa Float32

            # inferred type
            mat = LinearIsotropicElastic{3}(Float32(1.0), Float32(0.3))
            @test mat isa LinearIsotropicElastic{3, Float32}

            # mixed type
            mat = LinearIsotropicElastic{3}(1, 0.3) # Int + Float
            @test mat isa LinearIsotropicElastic{3, Float64}
        end
        
        # invalid parameters
        @test_throws ArgumentError LinearIsotropicElastic{3}(1.0, 0.5)   # nu = 0.5
        @test_throws ArgumentError LinearIsotropicElastic{3}(0.0, 0.3)     # E = 0
        @test_throws ArgumentError LinearIsotropicElastic{2}(1.0, -1.0)  # nu = -1
        @test_throws ArgumentError LinearIsotropicElastic{4}(1.0, 0.3)   # invalid dim
    end

    @testset "lame_constants" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.3)
        μ = mat.E / (2 * (1 + mat.nu))
        λ = mat.E * mat.nu / ((1 + mat.nu) * (1 - 2mat.nu))
        @test FFTMechHomo.lame_constants(mat.E, mat.nu) == (μ, λ)
    end

    @testset "compute_stress!" begin
        @testset "2D" begin
            mat = LinearIsotropicElastic{2}(1.0, 0.3)
            μ, λ = FFTMechHomo.lame_constants(mat.E, mat.nu)

            # ε₁₁ ≠ 0
            strain = zeros(3, 1, 1)
            stress = zeros(3, 1, 1)
            strain[1, 1, 1] = 1.0
            FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1))

            @test stress[1, 1, 1] ≈ (2μ + λ) * strain[1, 1, 1] # σ₁₁
            @test stress[2, 1, 1] ≈ λ * strain[1, 1, 1] # σ₂₂
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
            mat = LinearIsotropicElastic{3}(1.0, 0.3)
            μ, λ = FFTMechHomo.lame_constants(mat.E, mat.nu)

            # ε₁₁ ≠ 0
            strain = zeros(6, 1, 1, 1)
            stress = zeros(6, 1, 1, 1)
            strain[1, 1, 1, 1] = 1.0
            FFTMechHomo.compute_stress!(stress, strain, mat, CartesianIndex(1, 1, 1))

            @test stress[1, 1, 1, 1] ≈ (2μ + λ) * strain[1, 1, 1, 1] # σ₁₁
            @test stress[2, 1, 1, 1] ≈ λ * strain[1, 1, 1, 1] # σ₂₂
            @test stress[3, 1, 1, 1] ≈ λ * strain[1, 1, 1, 1] # σ₃₃
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