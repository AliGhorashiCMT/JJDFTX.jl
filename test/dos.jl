@testset "Make Mesh" begin
    mesh_array = Vector{Float64}[]
    for (i, j) in Tuple.(CartesianIndices(rand(10, 10)))
        push!(mesh_array, [(j-1)/10, (i-1)/10, 0])
    end
    @test permutedims(hcat(mesh_array...), (2, 1))  == JJDFTX.make_mesh(10, Val(2))
    mesh_array = Vector{Float64}[]
    for (i, j, k) in Tuple.(CartesianIndices(rand(10, 10, 10)))
        push!(mesh_array, [(k-1)/10, (j-1)/10, (i-1)/10])
    end
    permutedims(hcat(mesh_array...), (2, 1)) == JJDFTX.make_mesh(10, Val(3))
end


@testset "Alumnium" begin
    dir = "../data/electron_phonon_reference/"
    @test JJDFTX.load_bands_points("$dir/bandstruct.eigenvals", "$dir/bandstruct.kpoints") == (69, 10)
    @test length(bandstructkpoints2q(kpointsfile="$dir/bandstruct.kpoints")) == 69
    export_hwannier("$dir/wannier", [12, 12, 12], spin = Val('n'))
    Hwannier, cell_map = hwannier("$dir/wannier"), np.loadtxt("$dir/wannier.map.txt")
    E1, D1 = density_of_states(Hwannier, cell_map, Val(3), monte_carlo=true, histogram_width=10, 
    mesh=10, num_blocks=250, degeneracy=2)
    E2, D2 = density_of_states("$dir/totalE.dos")
    indx1 = argmin(abs.(E1 .- 4))
    indx2 = argmin(abs.(E2 .- 4))

    @test (D1[indx1] - D2[indx2])/D2[indx2]*100 < 8
    @test isapprox(find_chemical_potential(E1, D1)[2][end], 10, atol=1e-1)

end

#=@testset "graphene dos normalization" begin
    ##Test that DOS data taken directly from JDFTX has the correct (per unit cell) normalization. 
    ##Note that there are 8 bands, therefore the total integrated density of states should give us 16, when accounting for spin degeneracy.
    @test isapprox(graphene_dos_check(), 16, atol=1e-2) 
    @test isapprox(JJDFTX.check_graphene_dos(2.8, 1000, 100), 2, atol=1e-2) ## Tight binding description 
    
    Hwannier, cell_map = hwannier("wannierDefectBN33BCUp"),  np.loadtxt("wannierDefectBN33BCUp.map.txt")
    energies, dos = density_of_states(Hwannier, cell_map, num_blocks=10, mesh=36, histogram_width=1000, )

    xenergies, yoccupations = find_chemical_potential(HWannierdefect, cellmapdefect, mesh=50, histogram_width=1000, energy_range=1, offset=3, plotoccupations=false)
    onequarter  = xenergies[argmin(abs.(yoccupations .- 0.25))]
    halffilling = xenergies[argmin(abs.(yoccupations .- 0.5))]
    threequarter = xenergies[argmin(abs.(yoccupations .- 0.75))]

    @test isapprox(onequarter, -2.535, atol=1e-1)

end

@testset "Check Histogramming Compatibility with Numerical Integration" begin
    dir = joinpath(@__DIR__, "../data/one_band_models/")
    export_hwannier(dir*"wannier", [20, 20, 1], spin=Val('n'))
    HWannier, cellmap = hwannier(dir*"wannier.txt", dir*"wannier.map.txt", 1), np.loadtxt( dir*"wannier.map.txt")
    @test isapprox(1, density_of_states_wannier_quad_check(HWannier, cellmap, -3, 0, 250, δ=0.05, maxevals=2000), atol=5e-2)
end

@testset "Phonon Density of States" begin
    dir = joinpath(@__DIR__, "../data/Sodium/")
    forcematrix, cellmap = phonon_force_matrix(dir*"Na-lattice8")
    phonon_density_of_states(forcematrix, cellmap; mesh=5, histogram_width = 100, energy_range = 2)
    @test isapprox(sum(phonon_density_of_states(forcematrix, cellmap, Val(3), mesh=3, histogram_width=1000))/1000, 3, atol=1e-3)
    forcematrix, cellmap = phonon_force_matrix(dir*"Na-lattice8.phononCellMap", dir*"Na-lattice8.phononOmegaSq" )
end

=#