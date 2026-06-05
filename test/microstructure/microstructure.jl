@testset "Microstructure" begin
    @testset "Constructor" begin
        @testset "2D" begin
            n = Tuple(rand(1:4, 2))

            # Fill with different parameters
            materials = Array{FFTMechHomo.AbstractMaterial}(undef, n...)
            for i in eachindex(materials)
                E = rand()
                nu = rand() * 1.5 - 1.0
                materials[i] = LinearIsotropicElastic{2}(E, nu)
            end
            microstructure = Microstructure(materials)

            @test ndims(microstructure) == 2
            @test size(microstructure) == n
        end

        @testset "3D" begin
            n = Tuple(rand(1:4, 3))

            # Fill with different parameters
            materials = Array{FFTMechHomo.AbstractMaterial}(undef, n...)
            for i in eachindex(materials)
                E = rand()
                nu = rand() * 1.5 - 1.0
                materials[i] = LinearIsotropicElastic{3}(E, nu)
            end
            microstructure = Microstructure(materials)

            @test ndims(microstructure) == 3
            @test size(microstructure) == n
        end
        
        @testset "Errors" begin
            n = (4, 4, 4)
            materials = Array{FFTMechHomo.AbstractMaterial}(undef, n...)

            # Mixed material dimensions
            materials = Array{FFTMechHomo.AbstractMaterial}(undef, n...)
            fill!(materials, LinearIsotropicElastic{3}(1.0, 0.3))
            materials[1, 1, 1] = LinearIsotropicElastic{2}(1.0, 0.3)
            @test_throws ArgumentError Microstructure(materials)

            # Mixed types
            materials = Array{FFTMechHomo.AbstractMaterial}(undef, n...)
            fill!(materials, LinearIsotropicElastic{3}(1.0, 0.3))
            materials[1, 1, 1] = LinearIsotropicElastic{3}(Float32(1.0), Float32(0.3))
            @test_throws ArgumentError Microstructure(materials)

            # Material array dimensions don't match material dimensions
            mat = LinearIsotropicElastic{2}(1.0, 0.3)
            @test_throws ArgumentError Microstructure(fill(mat, n...))
        end  
    end
end