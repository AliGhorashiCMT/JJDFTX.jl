@testset "graphene dos normalization" begin
    ##Test that DOS data taken directly from JDFTX has the correct (per unit cell) normalization. 
    ##Note that there are 8 bands, therefore the total integrated density of states should give us 16, when accounting for spin degeneracy.
    @test isapprox(graphene_dos_check(), 16, atol=1e-2) 
    @test isapprox(JJDFTX.check_graphene_dos(2.8, 1000, 100), 2, atol=1e-2) ## Tight binding description 
end