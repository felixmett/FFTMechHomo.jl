@testset "MacroscopicStrain" begin
    @testset "Constructor" begin
        # implicit dim
        # 3D from length 6
        bc = MacroscopicStrain([1e-3, 0.0, 0.0, 0.0, 0.0, 0.0])
        @test bc isa MacroscopicStrain{3, Float64}
        @test bc.data == [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0]

        # 2D from length 3
        bc2 = MacroscopicStrain([1e-3, 0.0, 0.0])
        @test bc2 isa MacroscopicStrain{2, Float64}

        # Float32
        bc32 = MacroscopicStrain(Float32[1e-3, 0.0, 0.0, 0.0, 0.0, 0.0])
        @test bc32 isa MacroscopicStrain{3, Float32}

        # explicit dim
        bc = MacroscopicStrain{3}([1e-3, 0.0, 0.0, 0.0, 0.0, 0.0])
        @test bc isa MacroscopicStrain{3, Float64}

        bc2 = MacroscopicStrain{2}([1e-3, 0.0, 0.0])
        @test bc2 isa MacroscopicStrain{2, Float64}

        @testset "Errors" begin
            # invalid dim
            @test_throws ArgumentError MacroscopicStrain{4}([1e-3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

            # wrong number of components for dim
            @test_throws ArgumentError MacroscopicStrain{3}([1e-3, 0.0, 0.0])
            @test_throws ArgumentError MacroscopicStrain{2}([1e-3, 0.0, 0.0, 0.0, 0.0, 0.0])

            # invalid length for inferred dim
            @test_throws ArgumentError MacroscopicStrain([1e-3, 0.0, 0.0, 0.0])
        end
    end

    @testset "initialize_strain_field" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.3)
        ms = Microstructure(fill(mat, 4, 4, 4))
        ims = FFTMechHomo.InternalMicrostructure(ms)
        bc = MacroscopicStrain([1e-3, 0.0, 0.0, 0.0, 0.0, 0.0])

        strain = FFTMechHomo.initialize_strain_field(bc, ims)

        @test size(strain) == (6, 4, 4, 4)
        @test all(strain[1, :, :, :] .≈ 1e-3)
        @test all(strain[2:end, :, :, :] .≈ 0.0)
    end
end