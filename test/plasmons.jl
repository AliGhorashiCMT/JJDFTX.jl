
#=
We check the validity of plasmon dispersion methods 

=#
@testset "Intraband Plasmon Dispersion in Graphene" begin
    dir = "../data/RPA_Dielectric/"
    HWannier, cellmap = hwannier(dir*"wannier.txt", dir*"wannier.map.txt", 1), np.loadtxt(dir*"wannier.map.txt")
    plasmon = zeros(50)
    plasmon2 = zeros(50)
    lat = loadlattice(dir*"graphene.out")[2]
    for j in 1:50
        println(j)
        plasmon[j]=direct_epsilon_cubature(HWannier, cellmap, lat, [4/6*10/20, 0, 0], 8*j/50, -2, maxevals=6000)
    end
    plot(plasmon)

    impol = im_polarization(HWannier, cellmap, lat, [4/6*10/20, 0, 0],  -2, spin=1, mesh=200, histogram_width=20, normalized=false) 
    for j in 1:50
        println(j)
        plasmon2[j]=real(return_2d_epsilon([4/6*10/20, 0, 0], lat, 8*j/50, impol, 100, 20, false))
    end
    @test abs(argmin(log.(abs.(plasmon2))) - argmin(log.(abs.(plasmon))))  < 2

end