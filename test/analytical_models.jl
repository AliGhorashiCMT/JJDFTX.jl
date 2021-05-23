@testset "graphene conductivity" begin
    @test isapprox(1, imag(JJDFTX.graphene_conductivity(1, 1/50, 3)), atol=1e-1)
    @test isapprox(0, imag(JJDFTX.graphene_conductivity(1, 1/10, 1, atol=1e-10, delta=0.001)), atol=1e-2)
    @test  JJDFTX.exact_graphene_plasmon(1/6, 1) >  JJDFTX.exact_graphene_plasmon(1/12, 1)
end