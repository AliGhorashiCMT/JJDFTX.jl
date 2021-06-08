
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

    impol = im_polarization(HWannier, cellmap, lat, [4/6*10/20, 0, 0],  -2, spin=1, mesh=200, histogram_width=20, normalized=false) 
    for j in 1:50
        println(j)
        plasmon2[j]=real(return_2d_epsilon([4/6*10/20, 0, 0], lat, 8*j/50, impol, 100, 20, false))
    end
    @test abs(argmin(log.(abs.(plasmon2))) - argmin(log.(abs.(plasmon))))  < 2

end

#=
We check the plasmon dispersion of graphene with a one band model. 
=#
@testset "Graphene Drude Plasmon" begin
    plasmon = zeros(20, 20)

    dir = "../data/RPA_Dielectric/"
    lat = loadlattice(dir*"graphene.out")[2]

    HWannier, cellmap = hwannier(dir*"wannier.txt", dir*"wannier.map.txt", 1), np.loadtxt(dir*"wannier.map.txt")

    impols = []
    for (i, q) in enumerate(range(1/600, 1/6, length=20))
        #Note that for our graphene model, we must include the spin degeneracy
        impol = im_polarization(HWannier, cellmap, 1, 0, lat, [q, 0, 0],  -3.2, spin=2, mesh=200, histogram_width=20, normalized=false) 
        push!(impols, impol)
    end
    for (i, q) in enumerate(range(1/600, 1/6, length=20))
        println(i)
        impol = impols[i]
        for (j, ω) in enumerate(range(0.1, 3, length=20))
            plasmon[i, j]=real(return_2d_epsilon_scipy(q, ω, impol, 100, 20, 20))
        end
    end
    @test 1.2 < range(0.1, 3, length=20)[argmin(log.(abs.(plasmon[10, :])))] < 1.8
end

@testset "Graphene Drude Plasmon With ImPolarizationAtFilling" begin
    plasmon = zeros(20, 20)

    dir = "../data/RPA_Dielectric/"
    lat = loadlattice(dir*"graphene.out")[2]
    kpointsfile = dir*"bandstruct.kpoints"
    kpointsfile2 = dir*"bandstruct2.kpoints"
    HWannier, cellmap = hwannier(dir*"wannier.txt", dir*"wannier.map.txt", 1), np.loadtxt(dir*"wannier.map.txt")

    #We are going to do calculations at Fermi energy of 1eV with respect to the Dirac point. We first calculate what filling this corresponds to:
    #Note that we take into account the valley degeneracy but not the spin degeneracy
    filling = 2*π*(1/6)^2
    filling *= unit_cell_area(lat)/(4*π^2)
    #furthermore, each reciprocal lattice vector is about 3 inverse angstroms in length. The fermi wavevector is 1/6 inverse angstroms. 
    #So we need to go on a kpath from 0 to 0.06 of the length of one of the reciprocal lattice vectors
    #We use a separate bandstruct2.kpoints file to read in the relevant kpoints and interpolate them adequately.
    impols = im_polarizationatfilling(HWannier, cellmap, lat, filling, mesh=250, offset=5, histogram_width=10, energy_range=10, kpointsfile=kpointsfile2, spin=2)
    pl = return_2d_epsilons(0:0.01:3, impols, lat, max_energy=100, histogram_width=10, kpointsfile=kpointsfile2)
    @test 2 < (0:0.01:3)[argmin(log.(abs.(smooth(pl[61, :]))))] < 2.5 #Compare with known graphene plasmon relation using only intraband (Drude contribution)
end

#=
@testset "Wannier Landau Damping" begin
    dir = "../data/RPA_Dielectric/"
    lat = loadlattice(dir*"graphene.out")[2]
    kpointsfile = dir*"bandstruct.kpoints"
    kpointsfile2 = dir*"bandstruct2.kpoints"
    HWannier, cellmap = hwannier(dir*"wannier.txt", dir*"wannier.map.txt", 1), np.loadtxt(dir*"wannier.map.txt")
    losses = Float64[]
    for q in (0.1/6):(0.1/6):2/6
        #q = 2π*2.5*137/(4*pi*ħ*JJDFTX.c*1)*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units
        ω = JJDFTX.exact_graphene_plasmon(q, 1)
        println(q)
        loss=landau_damping(HWannier, cellmap, lat, 10, 60, [q, 0, 0], -3.3, 10)[round(Int, 10*ω)]
        println(loss)
        push!(losses, loss*ħ)
    end
end

=#