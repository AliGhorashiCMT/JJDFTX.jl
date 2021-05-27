@testset "graphene conductivity" begin
    @test isapprox(1, imag(JJDFTX.graphene_conductivity(1, 1/50, 3)), atol=1e-1)
    @test isapprox(0, imag(JJDFTX.graphene_conductivity(1, 1/10, 1, atol=1e-10, delta=0.001)), atol=1e-2)
    @test  JJDFTX.exact_graphene_plasmon(1/6, 1) >  JJDFTX.exact_graphene_plasmon(1/12, 1)

    #Below, we will check Kramers Kronig Relations for the imaginary and real polarizations at 0 doping for graphene. 
    q = 1/6
    a = kramers_kronig(ω->JJDFTX.imag_neutral(q, ω), 0.2,  60000, min_energy_integration=1.0000000000001) 
    b = JJDFTX.real_neutral(q, 0.2)
    @test abs((a-b)/a*100) < 1 #Check for less than 1 percent difference

    a = kramers_kronig_reverse(ω->JJDFTX.real_neutral(q, ω), 3, .99999999999)
    b = JJDFTX.imag_neutral(q, 3)
    @test abs((a-b)/a*100)  < 1

    q = 0.05
    a = kramers_kronig(ω->JJDFTX.imag_neutral(q, ω), .06,  6000000, min_energy_integration=.3000001)
    b = JJDFTX.real_neutral(q, 0.06)
    @test abs((a-b)/a*100) < 1 #Check for less than 1 percent difference

end