
@testset "Lattice Heat Capacity " begin
    #Below we verify the methods that are used to calculate the heat capactiy contribution from phonon 
    dir = "../data/electron_phonon_matrix_elements/"
    allattice = [[2, 2, 0], [0, 2, 2], [2, 0, 2]]
    forcematrix, cellmapph = phonon_force_matrix(joinpath(dir, "totalE"))
    wannierfile = joinpath(dir, "wannier.txt")
    cellmapfile = joinpath(dir, "wannier.map.txt")
    HWannier, cellmap = hwannier(wannierfile, cellmapfile, 5), np.loadtxt(cellmapfile)
    hcaps = lattice_heatcapacity.(10:100:2000, Ref(allattice), Ref(forcematrix), Ref(cellmapph), mesh=10)
    hcaps *= 1.60218e-19
    hcaps *= 1e30
    @test 20 <  hcaps[end]/(1e5) < 30 #Test that it comports with figure in Shankar's paper
end

@testset "Electronic Heat Capacity" begin
    #Below we verify the methods that are used to calculate the heat capacity contribution from electrons in Aluminum
    dir = "../data/electron_phonon_matrix_elements/"
    allattice = [[2, 2, 0], [0, 2, 2], [2, 0, 2]]
    μ = 10.88456
    forcematrix, cellmapph = phonon_force_matrix(joinpath(dir, "totalE"))
    wannierfile = joinpath(dir, "wannier.txt")
    cellmapfile = joinpath(dir, "wannier.map.txt")
    HWannier, cellmap = hwannier(wannierfile, cellmapfile, 5), np.loadtxt(cellmapfile)

    ecaps = electron_heatcapacity.(μ, range(500, 8000, length=10), Ref(allattice), Ref(HWannier), Ref(cellmap), 5, mesh=10, energy_range=40, histogram_width=5, offset=3)
    @test 7.4 < ecaps[end]*1.6*(1e-19)*(1e30)*(1e-5) < 8.1 #Check that it comports with figure in Shankar's paper
end 