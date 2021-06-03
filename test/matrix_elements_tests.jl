@testset "Electron-Phonon Matrix Element Tests" begin
    allattice = [[2, 2, 0], [0, 2, 2], [2, 0, 2]]
    μ = 0.4/eV
    #We test the equivalence of different methods of obtaining the Eliashbert spectral function and the methods therein
    dir = "../data/momentum_matrix_elements"
    wannierfile = joinpath(dir, "Al_wannierbands.txt")
    cellmapfile = joinpath(dir, "Al_cellmap.txt")
    HWannier, cellmap = hwannier(wannierfile, cellmapfile, 5), np.loadtxt(cellmapfile)
    dos = JJDFTX.density_of_states_montecarlo_3d(HWannier, cellmap, 5, mesh=35, energy_range=40, offset=10, histogram_width=40)
    dosmanuallyatmu = dos[round(Int, (μ+10)*40)]/unit_cell_volume(allattice)
    dosatmufunc = JJDFTX.dosatmu(HWannier, cellmap, allattice, 5, μ, mesh=35)
    dosatmulorentzian = JJDFTX.dosatmulorentzian(HWannier, cellmap, allattice, 5, μ, mesh=20, esmearing=0.1)
    dosatmugaussian = JJDFTX.dosatmugaussian(HWannier, cellmap, allattice, 5, μ, mesh=20, esmearing=0.1)
    @test abs(100*(dosmanuallyatmu-dosatmufunc)/dosatmufunc) < 5 #Less than five percent difference
    @test abs(100*(dosatmulorentzian-dosmanuallyatmu)/dosmanuallyatmu) < 5 #Less than five percent 
    @test abs(100*(dosatmugaussian-dosmanuallyatmu)/dosmanuallyatmu) < 5 #Less than five percent 
e