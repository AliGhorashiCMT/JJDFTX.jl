function dft_graphene_dos_per_area(;kwargs...)
    a = 1.42*sqrt(3)*1/bohrtoangstrom
    graphene_lattice = lattice([a -a/2 0; 0 a*sqrt(3)/2 0; 0 0 20])
    graphene_ucell_area = unit_cell_area(graphene_lattice)
    DOS_DATA_PATH = joinpath(@__DIR__, "../../data/graphene_examples/graphene.in.dos")
    plot(np.loadtxt(DOS_DATA_PATH)[:, 1]*27.2, np.loadtxt(DOS_DATA_PATH)[:, 2]/27.2/graphene_ucell_area, linewidth=4, 
        size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
end

function dft_graphene_phonon_dispersion(q::Vector{<:Real})
    dir =  "../../data/graphene_examples/"
    cell_map_path = joinpath(@__DIR__, dir*"graphene.in.phononCellMap")
    phonon_omegasq_path = joinpath(@__DIR__, dir*"graphene.in.phononOmegaSq")
    graphene_force_matrix, graphene_cell_map = phonon_force_matrix(cell_map_path, phonon_omegasq_path)
    return phonon_dispersion(graphene_force_matrix, graphene_cell_map, q)
end

"""
Returns the graphene energy dispersion along the Γ-M-K-Γ path in the Brillouin zone.
"""
function dft_graphene_wannier_dispersion()
    dir = "../../data/graphene_examples/"
    bands = zeros(8, 30)
    bands_dir = joinpath(@__DIR__, dir*"wannierbands.txt")
    map_dir = joinpath(@__DIR__, dir*"wanniercellmap.txt")
    for i in 1:10
        println(i);bands[:, i] = wannier_bands(bands_dir, map_dir, [0, 0.5*i/10, 0], 8)
    end
    for i in 1:10
        println(i+10);bands[:, i+10] = wannier_bands(bands_dir, map_dir, [2/3*i/10, 0.5-(0.5+1/3)*i/10, 0], 8)
    end
    for i in 1:10
        println(i+20);bands[:, i+20] = wannier_bands(bands_dir, map_dir, [2/3-2/3*i/10, -1/3+1/3*i/10, 0], 8)
    end
    return transpose(bands)
end


function dft_graphene_wannierbandsoverlayedDOS(mesh1::Integer=100, mesh2::Integer=1000)
    dir = "../../data/graphene_examples/"
    bands_dir = joinpath(@__DIR__, dir*"wannierbands.txt")
    map_dir = joinpath(@__DIR__, dir*"wanniercellmap.txt")
    DOS_DATA_PATH = joinpath(@__DIR__, dir*"graphene.in.dos")
    Hwannier = hwannier(bands_dir, map_dir, 8)
    cellmap = np.loadtxt(map_dir)
    bands = zeros(8, 300)
    for i in 1:100
        bands[:, i] = wannier_bands(Hwannier, cellmap, [0, 0.5*i/100, 0], 8)
    end
    for i in 1:100
        bands[:, i+100] = wannier_bands(Hwannier, cellmap,  [2/3*i/100, 0.5-(0.5+1/3)*i/100, 0], 8)
    end
    for i in 1:100
        bands[:, i+200] = wannier_bands(Hwannier, cellmap, [2/3-2/3*i/100, -1/3+1/3*i/100, 0], 8)
    end
    dosW = density_of_states_wannier(Hwannier, cellmap, 8, histogram_width = 5, mesh = mesh1, offset = 25, energy_range=35)

    C = plot(2*dosW[2], dosW[1], legend = false, title = "Wannier DOS", linewidth = 4)
    analyticgraphenedos = 2*graphene_dos(-2.8, mesh2, 5)
    D = plot(analyticgraphenedos,1:length(analyticgraphenedos), yticks = false, linewidth=4, title = "Dirac Cone Model DOS", legend = false)
    A = plot(transpose(bands), ylabel = "Energy (eV)", legend=false, linewidth = 4, xticks = false)
    B = plot(np.loadtxt(DOS_DATA_PATH)[:, 2]/27.2, np.loadtxt(DOS_DATA_PATH)[:, 1]*27.2, linewidth=4, yticks = false, legend = false, title="DFT DOS")
    plot(A, B, C,  D, size = (1000, 500))
    #return bands
end

function graphene_eph_matrix_elements(k1::Vector{<:Real}, k2::Vector{<:Real}, PhononBand::Integer)

    dir = "../../data/graphene_examples"
    bands_dir = joinpath(@__DIR__, dir*"/wannierbands.txt")
    map_dir = joinpath(@__DIR__, dir*"wanniercellmap.txt")

    cell_map_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellMap")
    cell_weights_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellWeights")

    cell_mapph_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellMapPh")
    cell_weightsph_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellWeightsPh")
    HePh_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfHePh")

    phonon_cellmap_dir = joinpath(@__DIR__, dir*"graphene.in.phononCellMap")
    phonon_omegasq_dir = joinpath(@__DIR__, dir*"graphene.in.PhononOmegaSq")

    HePhWannier, cellMapEph = write_eph_matrix_elements(cell_map_dir, cell_weights_dir, cell_mapph_dir, cell_weightsph_dir, HePh_dir, 6, [2, 2, 1])
    forcemat, mapph = phonon_force_matrix(phonon_cellmap_dir, phonon_omegasq_dir)

    return abs.(eph_matrix_elements(HePhWannier, cellMapEph, forcemat, mapph, bands_dir, map_dir, k1, k2, 8)[PhononBand, 4, :])
end

function graphene_eph_matrix_elements(HePhWannier, cellMapEph, forcemat, mapph, bands_dir, map_dir, k1, k2, PhononBand::Integer)
    return abs.(eph_matrix_elements(HePhWannier, cellMapEph, forcemat, mapph, bands_dir, map_dir, k1, k2, 8)[PhononBand, 4, :])
end



"""
$(TYPEDSIGNATURES)
We use an analytic model from the following: https://pubs.acs.org/doi/full/10.1021/nl402696q
In particular, we use the expressions given in table 1. Note that our comparisons here are in eV and we plot the square of the matrix element.

"""
function graphene_eph_matrix_elements_compare(numpoints::Integer)
    dir = "../../data/graphene_examples/"
    bands_dir = joinpath(@__DIR__, dir*"wannierbands.txt")
    map_dir = joinpath(@__DIR__, dir*"wanniercellmap.txt")
    cell_map_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellMap")
    cell_weights_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellWeights")
    cell_mapph_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellMapPh")
    cell_weightsph_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfCellWeightsPh")
    HePh_dir = joinpath(@__DIR__, dir*"wannier.graphene.in.mlwfHePh")
    phonon_cellmap_dir = joinpath(@__DIR__, dir*"graphene.in.phononCellMap")
    phonon_omegasq_dir = joinpath(@__DIR__, dir*"graphene.in.phononOmegaSq")
    HePhWannier, cellMapEph=write_eph_matrix_elements(cell_map_dir, cell_weights_dir, cell_mapph_dir, cell_weightsph_dir, HePh_dir, 6, [2, 2, 1])
    forcemat, mapph = phonon_force_matrix(phonon_cellmap_dir, phonon_omegasq_dir)
    K=[2/3, -1/3, 0]
    numerical5 = Vector{Float64}()
    numerical6 = Vector{Float64}()
    analytic5 = Vector{Float64}()
    analytic6 = Vector{Float64}()
    for i in 1:numpoints
        println("Progress: ", i)
        delta1 = [.05*i/(3*numpoints)+.03, .02*i/(3*numpoints)+.02, 0]
        delta2 = [.05-i*.001/(3*numpoints), .02-i*.02/(3*numpoints), 0]
        kvector1 = K + delta1
        kvector2 = K + delta2

        kvector=delta1 #Wave vector of initial electron with respect to K point
        qvector=delta2-delta1 ##Wave vector of phonon 
        kplusqvector=delta2 #Wave vector of final electron with respect to K point

        thetak = np.arctan(kvector[2]/kvector[1])
        thetaq = np.arctan(qvector[2]/qvector[1])

        thetakplusq = np.arctan(kplusqvector[2]/kplusqvector[1])
        theta1 = thetak-thetaq
        theta2 = thetakplusq-thetaq
            
        A1 = .045*(1+cos(theta1+theta2))
        A2 = .045*(1-cos(theta1+theta2))

        push!(numerical5, graphene_eph_matrix_elements(HePhWannier, cellMapEph, forcemat, mapph, bands_dir, map_dir, kvector1, kvector2, 5)[5]^2)
        push!(numerical6, graphene_eph_matrix_elements(HePhWannier, cellMapEph, forcemat, mapph, bands_dir, map_dir, kvector1, kvector2, 6)[5]^2)

        push!(analytic5, A1)
        push!(analytic6, A2)

    end
    plot(numerical5, label = "Numerical 5th Phonon Band", color="black", linewidth = 5)
    plot(numerical6, label = "Numerical 6th Phonon Band", color="black", linewidth = 5)
    plot(analytic5, label = "Analytic 5th Phonon Band", color="red", linewidth = 5)
    plot(analytic6, label = "Analytic 6th Phonon Band", color="red", linewidth = 5)
    ylabel("g^2(eV^2)")
    xticks(Float64[])

end

"""
A simple example to show the normalization of JDFTX data- This function should give a value close to 16 (8 bands and each has two possible spins)
"""
function graphene_dos_check()
    dir = "../../data/graphene_examples/"
    DOS_DATA_PATH = joinpath(@__DIR__, dir*"graphene.in.dos")
    x, y = np.loadtxt(DOS_DATA_PATH)[:, 1]*27.2, np.loadtxt(DOS_DATA_PATH)[:, 2]/27.2
    sum(y[2:end].*diff(x))
end

function graphene_wannier_impolarization(qx::Real; mesh::Integer = 20, histogram_width::Real = 10, degeneracy::Integer=4, subset::Integer=10, μ=-3)
    a = 1.42*sqrt(3)
    dir = "../../data/graphene_examples/"
    bands_dir = joinpath(@__DIR__, dir*"wannierbands.txt")
    map_dir = joinpath(@__DIR__, dir*"wanniercellmap.txt")
    HWannier=hwannier(bands_dir, map_dir, 8);
    cell_map=np.loadtxt(map_dir);
    im_pols = im_polarization(HWannier, cell_map, 8, 4, [[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]], [qx, 0, 0], μ; 
        spin=degeneracy, Koffset=[2/3, -1/3, 0], subset=subset, mesh=mesh, histogram_width=histogram_width) 
    return im_pols
end

"""
Given as an example of the general procedure to find the plasmon dispersion. First, the imaginary value of the polarization is found
"""
function example_graphene_wannier_plasmon(nqs::Integer, nomegas::Integer; mesh::Integer=30)
    plasmon = zeros(nqs, nomegas)
    for i in 1:nqs
        println(i)
        flush(stdout)
        gimpol = graphene_wannier_impolarization(i/nqs*1/6, mesh=mesh, histogram_width=100)
        for j in 1:nomegas
            plasmon[i, j] = return_2d_epsilon_scipy(i/nqs*1/6, 2*j/nomegas, smooth(gimpol, win_len=10), 100, 100, 30)
        end
    end
    return plasmon
end

function read_al_wannier_bands(k::Vector{<:Real})
    wannier_bands_path = joinpath(@__DIR__, "../../data/momentum_matrix_elements/Al_wannierbands.txt")
    cell_map_path = joinpath(@__DIR__, "../../data/momentum_matrix_elements/Al_cellmap.txt")
    Al_cellmap = np.loadtxt(cell_map_path)
    Al_hwannier = hwannier(wannier_bands_path, cell_map_path, 5)
    return wannier_bands(Al_hwannier, Al_cellmap, k, 5)
end

function example_al_wannier_bands()
    Al_BANDS = zeros(5, 500)
    for i in 1:100
        kx, ky, kz = 0,  0.5*i/100, 0.5*i/100
        Al_BANDS[:, i] = read_al_wannier_bands([kx, ky, kz])
    end
    for i in 1:100 
        kx, ky, kz = 0.25*i/100,  0.5+0.25*i/100, 0.5
        Al_BANDS[:, i+100] = read_al_wannier_bands([kx, ky, kz])
    end
    for i in 1:100 
        kx, ky, kz = 0.25+0.25*i/100,  0.75-0.25*i/100, 0.5
        Al_BANDS[:, i+200] = read_al_wannier_bands([kx, ky, kz])
    end
    for i in 1:100 
        kx, ky, kz = 0.5-0.5*i/100,  0.5-0.5*i/100, 0.5-0.5*i/100
        Al_BANDS[:, i+300] = read_al_wannier_bands([kx, ky, kz])
    end
    for  i in 1:100
        kx, ky, kz = 0.375*i/100, 0.375*i/100, 0.375*i/100
        Al_BANDS[:, i+400] = read_al_wannier_bands([kx, ky, kz])
    end
    plot([Al_BANDS[i, :] for i in 1:5])
    return Al_BANDS
end    

function example_aluminum_imepsilon(;histogram_width::Integer=10, mesh::Integer=10)
    dir = "../../data/momentum_matrix_elements/"
    wannier_bands_path = joinpath(@__DIR__, dir*"Al_wannierbands.txt")
    cell_map_path = joinpath(@__DIR__, dir*"Al_cellmap.txt")
    Pwannier_path = joinpath(@__DIR__, dir*"AlP.txt")
    Al_cellmap = np.loadtxt(cell_map_path)
    Al_hwannier = hwannier(wannier_bands_path, cell_map_path, 5)
    AlPwannier = pwannier(Pwannier_path, cell_map_path, 5)
    a = 4.05
    lattice_vectors = [[a/2, a/2, 0], [0, a/2, a/2], [a/2, 0, a/2]]
    return im_epsilon_3d_mc(lattice_vectors, Al_hwannier, Al_cellmap, AlPwannier, 5, 11; mesh=mesh, spin=2, histogram_width= histogram_width )
end

function get_RPA_dir()
    dir = "../../data/RPA_Dielectric/"

    return joinpath(@__DIR__, dir)
end

function get_Boltzmann_dir()
    dir = "../../data/boltzmann/"

    return joinpath(@__DIR__, dir)
end