@testset "Tensors" begin
    @testset "voigt_weights" begin
        # correct diagonal values
        @test diag(FFTMechHomo.voigt_weights(Val(3), Float64, Val(:strain))) == [1.0, 1.0, 0.5]
        @test diag(FFTMechHomo.voigt_weights(Val(6), Float64, Val(:strain))) == [1.0, 1.0, 1.0, 0.5, 0.5, 0.5]
        @test diag(FFTMechHomo.voigt_weights(Val(3), Float64, Val(:stress))) == [1.0, 1.0, 2.0]
        @test diag(FFTMechHomo.voigt_weights(Val(6), Float64, Val(:stress))) == [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
        @test diag(FFTMechHomo.voigt_weights(Val(3), Float64, Val(:mixed))) == [1.0, 1.0, 1.0]
        @test diag(FFTMechHomo.voigt_weights(Val(6), Float64, Val(:mixed))) == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        # stress and strain weights multiply to ones
        W_strain = FFTMechHomo.voigt_weights(Val(3), Float64, Val(:strain))
        W_stress = FFTMechHomo.voigt_weights(Val(3), Float64, Val(:stress))
        @test diag(W_strain) .* diag(W_stress) == [1.0, 1.0, 1.0]

        # invalid size
        @test_throws ArgumentError FFTMechHomo.voigt_weights(Val(4), Float64, Val(:strain))
    end

    @testset "tensor_dot" begin
        # no shear -> weights don't matter
        u = [1.0, 0.0, 0.0]
        @test FFTMechHomo.tensor_dot(u, u, :stress) ≈ 1.0
        @test FFTMechHomo.tensor_dot(u, u, :strain) ≈ 1.0

        # pure shear -> weight is 2 for stress, 0.5 for strain
        v = [0.0, 0.0, 1.0]
        @test FFTMechHomo.tensor_dot(v, v, :stress) ≈ 2.0
        @test FFTMechHomo.tensor_dot(v, v, :strain) ≈ 0.5

        # stress-strain mixed -> identity weights
        @test FFTMechHomo.tensor_dot(v, v, :mixed) ≈ 1.0

        # self tensor dot product form
        σ = [1.0, 2.0, 3.0]
        ε = [0.1, 0.2, 0.3]
        @test FFTMechHomo.tensor_dot(σ, :stress) ≈ FFTMechHomo.tensor_dot(σ, σ, :stress)
        @test FFTMechHomo.tensor_dot(ε, :strain) ≈ FFTMechHomo.tensor_dot(ε, ε, :strain)

        # unequal lengths
        @test_throws ArgumentError FFTMechHomo.tensor_dot([1.0, 0.0], [1.0, 0.0, 0.0], :stress)

        # Float32
        u32 = Float32[1.0, 0.0, 0.0]
        @test FFTMechHomo.tensor_dot(u32, u32, :stress) isa Float32
    end

    @testset "inner_product" begin
        # uniform field -> inner product equals tensor_dot at one point
        U = ones(3, 4, 4)
        volume = Float64(prod(size(U)[2:end]))
        @test FFTMechHomo.inner_product(U, volume, :stress) ≈ FFTMechHomo.tensor_dot(ones(3), :stress)

        # zero field
        @test FFTMechHomo.inner_product(zeros(3, 4, 4), volume, :stress) ≈ 0.0

        # self inner product form
        U = rand(3, 4, 4)
        @test FFTMechHomo.inner_product(U, U, volume, :stress) ≈ FFTMechHomo.inner_product(U, volume, :stress)
    end
end