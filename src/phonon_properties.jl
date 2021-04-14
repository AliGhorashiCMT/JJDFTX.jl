"Plots the phonon band dispersion at the kpoints supplied"
function plot_phonons(cell_map::String, phononOmegaSq::String, kpoints::String)
    cellMapPh = np.loadtxt(cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    kpointsIn = np.loadtxt(kpoints, skiprows=2, usecols=(1,2,3))
    nKin = size(kpointsIn)[1]
    #--- Fourier transform from real to k space:
    forceMatrixTilde = np.tensordot(np.exp((2im*np.pi)*np.dot(kpointsIn,transpose(cellMapPh))), forceMatrixPh, axes=1)
    #--- Diagonalize:
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTilde)
    plot(title="Phonon Dispersion", titlefontsize=20, ytickfontsize=15, sqrt.(abs.(omegaSq))/eV, linewidth=2, color="orange", legend=false, size=(800, 1000), xticks=[])
end

"Give phonon dispersion at individual kpoints"
function phonon_dispersion(phonon_cell_map::String, phononOmegaSq::String, qnorm::Array{<:Real, 1}) 
    cellMapPh = np.loadtxt(phonon_cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    #--- Fourier transform from real to k space:
    forceMatrixTildeq = np.tensordot(np.exp((2im*np.pi)*np.dot(qnorm,transpose(cellMapPh))), forceMatrixPh, axes=1)
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTildeq)
    return sqrt.(abs.(omegaSq))/eV
end

function phonon_dispersion(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorm::Array{<:Real, 1})
    forceMatrixTildeq = np.tensordot(np.exp(2im*π*np.dot(qnorm, transpose(phonon_cell_map)  )), force_matrix, axes=1   )
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTildeq)
    return sqrt.(abs.(omegaSq))/eV
end

function phonon_dispersionmodes(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorm::Array{<:Real, 1})
    phase = np.exp((2im*np.pi)*np.tensordot(qnorm, transpose(phonon_cell_map), axes=1))
    ### Note that we must permute the indices of the force matrix by jdftx convention
    omegaSq, U = np.linalg.eigh(np.tensordot(phase, permutedims(force_matrix, (1, 3, 2)), axes=1))
    return sqrt.(abs.(omegaSq))/eV, U
end

function phonon_force_matrix(phonon_cell_map::String, phononOmegaSq::String)
    cellMapPh = np.loadtxt(phonon_cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    return forceMatrixPh, cellMapPh
end

function phonon_force_matrix(filebase::String)
    phonon_cell_map = "$(filebase).phononCellMap"
    phononOmegaSq = "$(filebase).phononOmegaSq"
    cellMapPh = np.loadtxt(phonon_cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    return forceMatrixPh, cellMapPh
end

#=
function plot_phononsjl(cell_map::String, phononOmegaSq::String, kpoints::String)
    cellMapPh=readdlm(cell_map)[2:end,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))

    eV=1/27.2
    kpointsIn = np.loadtxt(kpoints, skiprows=2, usecols=(1,2,3))
    nKin = size(kpointsIn)[1]
    #--- Fourier transform from real to k space:
    forceMatrixTilde = np.tensordot(np.exp((2im*np.pi)*np.dot(kpointsIn,transpose(cellMapPh))), forceMatrixPh, axes=1)
    #--- Diagonalize:
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTilde)
    plot(title="Phonon Dispersion", titlefontsize=20, ytickfontsize=15, sqrt.(abs.(omegaSq))/eV, linewidth=2, color="orange", legend=false, size=(800, 1000), xticks=[])

end

=#

"""
Returns the electron self energy to lowest order in the electron-phonon interaction. The expression used is from the paper:
Park, Cheol-Hwan, et al. "Velocity renormalization and carrier lifetime in graphene from the electron-phonon interaction." Physical review letters 99.8 (2007): 086804.
This is commonly referred to as the Migdal approximation. 
The graphene methods for self energy use the same approximation
"""
function migdal_approximation(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real; histogram_length::Real=100, mesh::Int=30) 
    qnormalized = normalize_kvector(lattice_vectors, q)
    self_energy = 0
    ϵi = wannier_bands(HWannier, cell_map, qnormalized)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            phonon_energies = phonon_dispersion(force_matrix, phonon_cell_map, [xmesh/mesh, ymesh/mesh, 0])
            phonon_mat_elements= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, qnormalized, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
            ϵf = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
            fermi = ϵf<μ ? 1 : 0
            for phonon in 1:length(phonon_energies)
                ωph = phonon_energies[phonon] 
                if abs((ϵi-ϵf-ωph)*histogram_length)<0.5
                    self_energy +=  π*abs(phonon_mat_elements[phonon])^2*(1-fermi)*histogram_length/mesh^2
                end
                if abs((ϵi-ϵf+ωph)*histogram_length)<0.5
                    self_energy +=  π*abs(phonon_mat_elements[phonon])^2*(fermi)*histogram_length/mesh^2
                end
            end
        end
    end
end