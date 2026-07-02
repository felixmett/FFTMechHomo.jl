microstructure_27_voxels(mat1, mat2) = Microstructure([i==1 ? mat1 : mat2 for i in 1:3, j in 1:3, k in 1:3])

@testset "Integration" begin
    @testset "HistoryIndependent Materials" begin
        # reference average stress for 27 voxel microstructure
        # testcase source: https://fft-workshop-22.sciencesconf.org/ (Tutorial sessions)
        reference_27_voxels = [
            79.3500, 26.2562, 26.2562, 0.0, 0.0, 0.0
        ]

        mat_soft = LinearIsotropicElastic{3}(2100., 0.3)
        mat_hard = LinearIsotropicElastic{3}(72000., 0.22)
        ms = microstructure_27_voxels(mat_soft, mat_hard)

        α₀ = mat_soft.μ + mat_hard.μ
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
    
    @testset "Custom HistoryIndependent Material" begin
        struct TestCustomElasticity{dim, T <: AbstractFloat} <: HistoryIndependent{dim, T}
            E::T
            nu::T
            μ::T
            λ::T
            function TestCustomElasticity{dim, T}(E::T, nu::T) where {dim, T <: AbstractFloat}
                μ = E / (2 * (1 + nu))
                λ = E * nu / ((1 + nu) * (1 - 2nu))
                new{dim, T}(E, nu, μ, λ)
            end
        end

        function FFTMechHomo.compute_stress!(
            stress::AbstractArray{T},
            strain::AbstractArray{T},
            mat::TestCustomElasticity{dim, T},
            i::CartesianIndex
        ) where {dim, T <: AbstractFloat}
            tr_strain = sum(strain[1:dim, i])
            stress[1:dim, i] .= 2mat.μ .* strain[1:dim, i] .+ mat.λ * tr_strain

            stress[dim+1:end, i] .= mat.μ .* strain[dim+1:end, i]
            return
        end

        E_soft = 2100.
        nu_soft = 0.3
        E_hard = 72000.
        nu_hard = 0.22 
        dim = 3

        # Custom material
        mat_soft_custom = TestCustomElasticity{dim, Float64}(E_soft, nu_soft)
        mat_hard_custom = TestCustomElasticity{dim, Float64}(E_hard, nu_hard)
        microstructure = microstructure_27_voxels(mat_soft_custom, mat_hard_custom)

        disc = MoulinetSuquetDiscretization(microstructure)
        macro_strain = MacroscopicStrain([0.01, 0, 0, 0, 0, 0])
        solver = BasicScheme(mat_soft_custom.μ + mat_hard_custom.μ, microstructure)
        sol = solve(microstructure, disc, macro_strain, solver)

        # Equivalent built-in material
        mat_soft_builtin = LinearIsotropicElastic{dim, Float64}(E_soft, nu_soft)
        mat_hard_builtin = LinearIsotropicElastic{dim, Float64}(E_hard, nu_hard)
        validation_microstructure = microstructure_27_voxels(mat_soft_builtin, mat_hard_builtin)

        validation_disc = MoulinetSuquetDiscretization(validation_microstructure)
        validation_solver = BasicScheme(mat_soft_builtin.μ + mat_hard_builtin.μ, validation_microstructure)
        validation_sol = solve(validation_microstructure, validation_disc, macro_strain, validation_solver)

        @test sol.stress_avg ≈ validation_sol.stress_avg
    end
end