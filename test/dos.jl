@testset "graphene dos normalization" begin
    ##Test that DOS data taken directly from JDFTX has the correct (per unit cell) normalization. 
    ##Note that there are 8 bands, therefore the total integrated density of states should give us 16, when accounting for spin degeneracy.
    @test isapprox(graphene_dos_check(), 16, atol=1e-2) 
    @test isapprox(JJDFTX.check_graphene_dos(2.8, 1000, 100), 2, atol=1e-2) ## Tight binding description 
    cellmapdefect, HWannierdefect = np.loadtxt("wannierDefectBN33BCUp.map.txt"), hwannier("wannierDefectBN33BCUp.txt", "wannierDefectBN33BCUp.map.txt")
    xenergies, yoccupations = find_chemical_potential(HWannierdefect, cellmapdefect, mesh=50, histogram_width=1000, energy_range=1, offset=3, plotoccupations=false)
    onequarter  = xenergies[argmin(abs.(yoccupations .- 0.25))]
    halffilling = xenergies[argmin(abs.(yoccupations .- 0.5))]
    threequarter = xenergies[argmin(abs.(yoccupations .- 0.75))]

    @test isapprox(onequarter, -2.535, atol=1e-1)

end

@testset "Check Histogramming Compatibility with Numerical Integration" begin
    dir = joinpath(@__DIR__, "../data/one_band_models/")
    write_map_write_h(dir*"wannier", [20, 20, 1], spin=Val('n'))
    HWannier, cellmap = hwannier(dir*"wannier.txt", dir*"wannier.map.txt", 1), np.loadtxt( dir*"wannier.map.txt")
    @test wannier_vectors(HWannier, cellmap, rand(3)) ≈ [1] #Check that one band model indeed has no nontrivial band eigenvectors
    @test wannier_vectors(dir*"wannier.txt", dir*"wannier.map.txt", rand(3), 1) ≈ [1] #Check that one band model indeed has no nontrivial band eigenvectors
    @test isapprox(1, density_of_states_wannier_quad_check(HWannier, cellmap, -3, 0, 250, δ=0.05, maxevals=2000), atol=5e-2)
    @test isapprox(1, density_of_states_wannier_quad_check(dir*"wannier.txt", dir*"wannier.map.txt",  -3, 0, 250, δ=0.05, maxevals=2),  atol=5e-2)
end