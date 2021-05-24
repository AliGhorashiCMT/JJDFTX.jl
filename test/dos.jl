@testset "graphene dos normalization" begin
    ##Test that DOS data taken directly from JDFTX has the correct (per unit cell) normalization. 
    ##Note that there are 8 bands, therefore the total integrated density of states should give us 16, when accounting for spin degeneracy.
    @test isapprox(graphene_dos_check(), 16, atol=1e-2) 
    @test isapprox(JJDFTX.check_graphene_dos(2.8, 1000, 100), 2, atol=1e-2) ## Tight binding description 
    cellmapdefect, HWannierdefect = np.loadtxt("wannierDefectBN33BCUp.map.txt"), hwannier("wannierDefectBN33BCUp.txt", "wannierDefectBN33BCUp.map.txt")
    xenergies, yenergies = find_chemical_potential(HWannierdefect, cellmapdefect, mesh=50, histogram_width=1000, energy_range=1, offset=3, plotoccupations=true)
    onequarter  = xenergies[argmin(abs.(yoccupations .- 0.25))]
    halffilling = xenergies[argmin(abs.(yoccupations .- 0.5))]
    threequarter = xenergies[argmin(abs.(yoccupations .- 0.75))]

    @test isapprox(onequarter, -2.535, atol=1e-1)

end