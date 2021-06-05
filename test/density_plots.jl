@testset "hBN density and fft box" begin
    dir = "../data/plot_densities_wfns/"
    V, Box, N = plot_density(dir*"graphene_density.n", dir*"graphene_density.out" , Val('z'), 2)
    @test Box == [54, 54, 168]
    @test isapprox(V, 4.651*4.651*15*sqrt(3)/2, atol=1e-2)
    @test isapprox(N, 8, atol=1e-2)
    V, Box, N = plot_density(dir*"graphene_density.n", dir*"graphene_density.out" , Val('x'), 2)
    @test Box == [54, 54, 168]
    @test isapprox(N, 8, atol=1e-2)
    @test isapprox(V, 4.651*4.651*15*sqrt(3)/2, atol=1e-2)
    V, Box, N = plot_density(dir*"graphene_density.n", dir*"graphene_density.out" , Val('y'), 2)
    @test Box == [54, 54, 168]
    @test isapprox(V, 4.651*4.651*15*sqrt(3)/2, atol=1e-2)
    @test isapprox(N, 8, atol=1e-2)
end


@testset "spin polarized density and fft box " begin
    dir = "../data/plot_densities_wfns/"
    plot_diffdensity(dir*"spin_density.n_dn", dir*"spin_density.n_up",  dir*"spin_density.out" , Val('x'), 2, 30, legend=false, xticks=false, yticks=false )
    plot_diffdensity(dir*"spin_density.n_dn", dir*"spin_density.n_up",  dir*"spin_density.out" , Val('y'), 2, 30, legend=false, xticks=false, yticks=false )
    V, Box, NUp, NDn = plot_diffdensity(dir*"spin_density.n_dn", dir*"spin_density.n_up",  dir*"spin_density.out" , Val('z'), 2, 30, legend=false, xticks=false, yticks=false )
    @test isapprox(1, abs(NDn-NUp), atol=1e-2) #Check that spin polarization is correct
    @test isapprox(33, NUp+NDn, atol=1e-2) 
end