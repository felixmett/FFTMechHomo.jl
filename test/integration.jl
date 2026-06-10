function microstructure_27_voxels()
    mat1 = LinearIsotropicElastic{3}(2100., 0.3)
    mat2 = LinearIsotropicElastic{3}(72000., 0.22)
    materials = [
        i == 1 ? mat1 : mat2 for i in 1:3, j in 1:3, k in 1:3
    ]
    return Microstructure(materials)
end

@testset "Integration" begin
    # reference average stress for 27 voxel microstructure
    # testcase source: https://fft-workshop-22.sciencesconf.org/ (Tutorial sessions)
    reference_27_voxels = [
        79.3500, 26.2562, 26.2562, 0.0, 0.0, 0.0
    ]

    ms = microstructure_27_voxels()
    α₀ = 0.5 * (2100 + 72000)
    macro_strain = MacroscopicStrain([0.01, 0, 0, 0, 0, 0])

    discretizations = [
        MoulinetSuquetDiscretization(ms)
    ]
    solvers = [
        BasicScheme(α₀, ms)
    ]
    for disc in discretizations
        for solver in solvers
            @testset "$(typeof(disc)) + $(typeof(solver))" begin
                sol = solve(ms, disc, macro_strain, solver)
                @test isapprox(sol.stress_avg, reference_27_voxels, rtol=1e-4)
            end
        end
    end
end