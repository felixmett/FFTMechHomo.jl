@testset "BasicScheme" begin
    @testset "Constructor" begin
        mat = LinearIsotropicElastic{3}(1.0, 0.3)
        ms = Microstructure(fill(mat, 4, 4, 4))

        solver = BasicScheme(2 * mat.μ, ms)
        @test solver isa BasicScheme{FFTMechHomo.Linear, Float64}
        @test solver.α₀ ≈ 2 * mat.μ
        @test solver.tol ≈ 1e-5
        @test solver.maxiter == 100
        @test solver.FFTW_flags == FFTW.MEASURE
        @test solver.FFTW_num_threads == 1

        # keyword arguments
        solver2 = BasicScheme(2 * mat.μ, ms; tol=1e-6, maxiter=200, FFTW_num_threads=4)
        @test solver2.tol ≈ 1e-6
        @test solver2.maxiter == 200
        @test solver2.FFTW_num_threads == 4

        # type promotion — Int α₀ promoted to T from microstructure
        solver3 = BasicScheme(2, ms)
        @test solver3 isa BasicScheme{FFTMechHomo.Linear, Float64}
        @test solver3.α₀ isa Float64

        # Float32 microstructure
        mat32 = LinearIsotropicElastic{3, Float32}(Float32(1), Float32(0.3))
        ms32 = Microstructure(fill(mat32, 4, 4, 4))
        solver32 = BasicScheme(Float32(2 * mat.μ), ms32)
        @test solver32 isa BasicScheme{FFTMechHomo.Linear, Float32}
    end
end