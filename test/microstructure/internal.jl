@testset "InternalMicrostructure" begin

    @testset "Constructor — single material" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.3)
        ms = Microstructure(fill(mat, 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)

        @test ndims(ims) == 3
        @test size(ims) == (4, 4, 4)
        @test length(ims.groups) == 1
        @test ims.groups[1] isa FFTMechHomo.MaterialGroup{3, LinearIsotropicElastic{3, Float64}}
        @test length(ims.groups[1].indices) == 4^3
    end

    @testset "Constructor — two materials" begin
        mat1 = LinearIsotropicElastic{3}(1.0, 0.3)
        mat2 = FFTMechHomo.ReferenceMaterial{3}(2mat1.μ)
        materials = Array{FFTMechHomo.AbstractMaterial}(undef, 4, 4, 4)
        fill!(materials, mat1)
        materials[1, 1, 1] = mat2

        ms = Microstructure(materials)
        ims = FFTMechHomo.InternalMicrostructure(ms)

        @test length(ims.groups) == 2
        total_indices = sum(length(g.indices) for g in ims.groups)
        @test total_indices == 4^3
    end

    @testset "compute_stress_field!" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.3)
        ms = Microstructure(fill(mat, 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)

        strain = zeros(6, 4, 4, 4)
        stress = zeros(6, 4, 4, 4)
        strain[1, :, :, :] .= 1e-3

        FFTMechHomo.compute_stress_field!(stress, strain, ims)

        # every voxel should have the same stress since material is uniform
        μ = mat.μ
        λ = mat.λ
        @test all(stress[1, :, :, :] .≈ (2μ + λ) * 1e-3)
        @test all(stress[2, :, :, :] .≈ λ * 1e-3)
        @test all(stress[3, :, :, :] .≈ λ * 1e-3)
        @test all(stress[4:6, :, :, :] .≈ 0.0)
    end

    @testset "compute_polarization_field!" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.0)
        ms = Microstructure(fill(mat, 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)
        ref = FFTMechHomo.ReferenceMaterial{3}(mat.μ * 2)

        strain = zeros(6, 4, 4, 4)
        stress = zeros(6, 4, 4, 4)
        strain[1, :, :, :] .= 1e-3

        FFTMechHomo.compute_polarization_field!(stress, strain, ims, ref)

        # polarization = σ(ε) - σ₀(ε) should be zero everywhere if the materials are equal
        @test all(stress .≈ 0)
    end

end