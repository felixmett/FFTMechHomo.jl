@testset "Solver utilities" begin
    @testset "GridConstants" begin
        ms = Microstructure(fill(LinearIsotropicElastic{3}(210e3, 0.3), 8, 8, 8))
        ims = FFTMechHomo.InternalMicrostructure(ms)
        consts = FFTMechHomo.GridConstants(ims)

        @test consts.n_voigt == 6
        @test consts.n_cells == 8^3
        @test consts.zero_idx == CartesianIndex(1, 1, 1)

        # 2D
        ms2 = Microstructure(fill(LinearIsotropicElastic{2}(210e3, 0.3), 8, 8))
        ims2 = FFTMechHomo.InternalMicrostructure(ms2)
        consts2 = FFTMechHomo.GridConstants(ims2)

        @test consts2.n_voigt == 3
        @test consts2.n_cells == 8^2
        @test consts2.zero_idx == CartesianIndex(1, 1)
    end

    @testset "cell_average!" begin
        ms = Microstructure(fill(LinearIsotropicElastic{3}(210e3, 0.3), 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)
        consts = FFTMechHomo.GridConstants(ims)

        # construct a field and its FFT
        field = ones(6, 4, 4, 4)
        field_fft = rfft(field, 2:4)
        cell_avg = zeros(6)

        FFTMechHomo.cell_average!(cell_avg, field_fft, consts)

        # average of a uniform field of ones should be ones
        @test all(cell_avg .≈ 1.0)
    end

    @testset "impose_macroscopic_strain!" begin
        ms = Microstructure(fill(LinearIsotropicElastic{3}(210e3, 0.3), 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)
        consts = FFTMechHomo.GridConstants(ims)

        strain_macro = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0]
        bc = FFTMechHomo.MacroscopicStrain(strain_macro)
        field_fft = zeros(Complex{Float64}, 6, 4, 4, 3)

        FFTMechHomo.impose_macroscopic_strain!(field_fft, bc, consts)

        # zero frequency mode should be strain_macro * n_cells
        @test all(real(field_fft[:, 1, 1, 1]) .≈ strain_macro * consts.n_cells)
    end

    @testset "residual!" begin
        mat = LinearIsotropicElastic{3}(210e3, 0.3)
        ms = Microstructure(fill(mat, 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)
        ref = FFTMechHomo.ReferenceMaterial{3}(mat.μ * 2)
        consts = FFTMechHomo.GridConstants(ims)

        strain = zeros(6, 4, 4, 4)
        strain[1, :, :, :] .= 1e-3
        strain_prev = zero(strain)  # previous strain is zero
        stress = zero(strain)
        stress_avg = ones(6)  # nonzero to avoid division by zero

        r = FFTMechHomo.residual!(stress, strain, strain_prev, stress_avg, ref, consts)

        @test r isa Float64
        @test r > 0.0
        @test isfinite(r)

        # zero strain difference should give zero residual
        r_zero = FFTMechHomo.residual!(stress, strain, strain, stress_avg, ref, consts)
        @test r_zero ≈ 0.0
    end
end