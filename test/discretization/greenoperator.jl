function test_discretizations(n::Int, dim::Int)
    implemented_discs = [
        MoulinetSuquetDiscretization(Microstructure(fill(LinearIsotropicElastic{dim}(1.0, 0.3), ntuple(x -> n, dim)...)))
    ]
    return implemented_discs
end

@testset "AbstractDiscreteGreenOperator" begin
    @testset "is_zero_or_nyquist" begin
        @testset "2D" begin
            for n in (7, 8)
                # only even grid lengths have nyquist frequencies
                has_nyquist = n % 2 == 0 ? 1 : 0
                for disc in test_discretizations(n, 2)
                    # zero frequency
                    @test FFTMechHomo.is_zero_or_nyquist(CartesianIndex(1, 1), disc)

                    # number of nyquist frequency indices and zero frequency should be flagged
                    n_flagged_target = has_nyquist * (2 * (n - 1) + 1) + 1
                    n_flagged = count(idx -> FFTMechHomo.is_zero_or_nyquist(idx, disc), CartesianIndices(disc.grid_size))
                    @test n_flagged == n_flagged_target
                end
            end
        end

        @testset "3D" begin
            for n in (7, 8)
                # only even grid lengths have nyquist frequencies
                has_nyquist = n % 2 == 0 ? 1 : 0
                for disc in test_discretizations(n, 3)
                    # zero frequency
                    @test FFTMechHomo.is_zero_or_nyquist(CartesianIndex(1, 1, 1), disc)

                    # number of nyquist frequency indices and zero frequency should be flagged
                    n_flagged_target = has_nyquist * (3 * (n - 1)^2 + 3 * (n - 1) + 1) + 1
                    n_flagged = count(idx -> FFTMechHomo.is_zero_or_nyquist(idx, disc), CartesianIndices(disc.grid_size))
                    @test n_flagged == n_flagged_target
                end
            end
        end
    end

    @testset "voigt_index" begin
        @testset "2D" begin
            voigt_idx2 = FFTMechHomo.voigt_index(Val(2))

            # diagonal entries
            @test voigt_idx2[1][1] == 1   # (1,1) -> 1
            @test voigt_idx2[2][2] == 2   # (2,2) -> 2

            # off-diagonal symmetry
            @test voigt_idx2[1][2] == voigt_idx2[2][1] == 3   # (1,2) -> 3, (2,1) -> 3

            # all Voigt indices covered
            @test sort(unique(collect(Iterators.flatten(voigt_idx2)))) == [1, 2, 3]
        end

        @testset "3D" begin
            voigt_idx3 = FFTMechHomo.voigt_index(Val(3))

            # diagonal entries
            @test voigt_idx3[1][1] == 1
            @test voigt_idx3[2][2] == 2
            @test voigt_idx3[3][3] == 3

            # off-diagonal symmetry
            @test voigt_idx3[2][3] == voigt_idx3[3][2]  # (2,3) -> 4, (3,2) -> 4
            @test voigt_idx3[1][3] == voigt_idx3[3][1]  # (1,3)-> 5, (3,1) -> 5
            @test voigt_idx3[1][2] == voigt_idx3[2][1]  # (1,2) -> 6, (2,1) -> 6

            # all Voigt indices covered
            @test sort(unique(collect(Iterators.flatten(voigt_idx3)))) == [1, 2, 3, 4, 5, 6]
        end
    end

    @testset "rfft_output_size" begin
        @testset "2D" begin
            for n in (7, 8)
                disc = MoulinetSuquetDiscretization(Microstructure(fill(LinearIsotropicElastic{2}(1.0, 0.3), n, n)))

                # first dimension is n_voigt
                @test FFTMechHomo.rfft_output_size(disc)[1] == 3

                # first spatial dimension is halved due to rfft
                expected_size = (3, n ÷ 2 + 1, n)
                rfft_output_size = FFTMechHomo.rfft_output_size(disc)
                @test rfft_output_size == rfft_output_size
            end
        end

        @testset "3D" begin
            for n in (7, 8)
                disc = MoulinetSuquetDiscretization(Microstructure(fill(LinearIsotropicElastic{3}(1.0, 0.3), n, n, n)))

                # first dimension is n_voigt
                @test FFTMechHomo.rfft_output_size(disc)[1] == 6

                # first spatial dimension is halved due to rfft
                expected_size = (6, n ÷ 2 + 1, n, n)
                rfft_output_size = FFTMechHomo.rfft_output_size(disc)
                @test rfft_output_size == rfft_output_size
            end
        end
    end
end