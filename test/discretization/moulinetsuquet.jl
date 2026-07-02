@testset "MoulinetSuquetDiscretization" begin
    @testset "Constructor" begin
        n = 8
        @testset "2D" begin
            ms = Microstructure(fill(LinearIsotropicElastic{2}(1.0, 0.3), n, n))
            disc = MoulinetSuquetDiscretization(ms)

            @test disc.grid_size == (n, n)
            @test ndims(disc) == 2
            @test eltype(disc) == Float64

            # frequency vectors
            @test length(disc.ξ) == 2
            @test length(disc.ξ[1]) == n ÷ 2 + 1
            @test length(disc.ξ[2]) == n
        end

        @testset "3D" begin
            ms = Microstructure(fill(LinearIsotropicElastic{3}(1.0, 0.3), n, n, n))
            disc = MoulinetSuquetDiscretization(ms)

            @test disc.grid_size == (n, n, n)
            @test ndims(disc) == 3
            @test eltype(disc) == Float64

            # frequency vectors
            @test length(disc.ξ) == 3
            @test length(disc.ξ[1]) == n ÷ 2 + 1
            @test length(disc.ξ[2]) == n
            @test length(disc.ξ[3]) == n
        end
    end
end