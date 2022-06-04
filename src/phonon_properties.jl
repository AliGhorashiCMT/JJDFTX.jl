"""
$(TYPEDSIGNATURES)

"""
function phonon_force_matrix(phonon_cell_map::AbstractString, phononOmegaSq::AbstractString)
    cellph_map = np.loadtxt(phonon_cell_map)[:,1:3]
    force_matrix = np.fromfile(phononOmegaSq, dtype=np.float64)
    num_cells = first(size(cellph_map))
    num_modes = Int(np.sqrt(first(size(force_matrix)) / num_cells))
    force_matrix = np.reshape(force_matrix, (num_cells, num_modes, num_modes))
    force_matrix = permutedims(force_matrix, (1, 3, 2))
    return force_matrix, cellph_map
end

phonon_force_matrix(filebase::AbstractString) = phonon_force_matrix("$(filebase).phononCellMap", "$(filebase).phononOmegaSq")

"""
$(TYPEDSIGNATURES)
"""
function diagonalize_phonon(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, k::Vector{<:Real})
    phase = np.exp(2im*π*cellph_map*k)
    omegaSq, U = np.linalg.eigh(np.tensordot(phase, force_matrix, axes=1))
    return omegaSq, U
end

function diagonalize_phonon(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, kpoints::AbstractArray{<:Real, 2})
    phase = np.exp(2im*π*cellph_map*kpoints)
    omegaSqs, Us = np.linalg.eigh(np.einsum("kij, kl -> lij", force_matrix, phase))
    return omegaSqs, Us
end

function diagonalize_phonon(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, k1s::AbstractArray{<:Real, 2}, 
    k2s::AbstractArray{<:Real, 2})
    kpoints = np.repeat(np.reshape(k1s, (3, -1, 1)), size(k2s)[2], axis=2) - 
            np.repeat(np.reshape(k2s, (3, 1, -1)), size(k1s)[2], axis=1)
    phase = np.exp(2im*π*np.tensordot(cellph_map, kpoints, axes=1))
    omegaSqs, Us = np.linalg.eigh(np.einsum("kij, klq -> lqij", force_matrix, phase))
    return omegaSqs, Us
end

"""
$(TYPEDSIGNATURES)

Optional keyword argument return_negative may be supplied to return imaginary frequencies as negative frequencies

"""
function phonon_dispersion(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, k::Vector{<:Real}; return_negative::Bool=false)
    omegaSq, _ = diagonalize_phonon(force_matrix, cellph_map, k)
    freq = return_negative ? sign.(omegaSq).*sqrt.(abs.(omegaSq))/eV : sqrt.(abs.(omegaSq))/eV
    return freq
end

function phonon_dispersion(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, kpoints::AbstractArray{<:Real, 2}; return_negative::Bool=false)
    omegaSq, _ = diagonalize_phonon(force_matrix, cellph_map, kpoints)
    freq = return_negative ? sign.(omegaSq).*sqrt.(abs.(omegaSq))/eV : sqrt.(abs.(omegaSq))/eV
    return freq
end

phonon_dispersion(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}; kpointsfile::AbstractString="bandstruct.kpoints", kwargs...) =
phonon_dispersion(force_matrix, cellph_map, hcat(bandstructkpoints2q(kpointsfile=kpointsfile)...); kwargs...)

"""
Returns the electron self energy to lowest order in the electron-phonon interaction. The expression used is from the paper:
Park, Cheol-Hwan, et al. "Velocity renormalization and carrier lifetime in graphene from the electron-phonon interaction." Physical review letters 99.8 (2007): 086804.
This is commonly referred to as the Migdal approximation. 
The graphene methods for self energy use the same approximation
"""
function migdal_approximation(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5},
    cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2},
    lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; histogram_width::Real=100, mesh::Integer=30) 

    qnormalized = normalize_kvector(lattice_vectors, q)
    self_energy = 0
    ϵi = wannier_bands(HWannier, cell_map, qnormalized)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        phonon_energies = phonon_dispersion(force_matrix, phonon_cell_map, [xmesh/mesh, ymesh/mesh, 0])
        phonon_mat_elements= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, qnormalized, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
        ϵf = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
        fermi = ϵf<μ ? 1 : 0
        for (idx, ωph) in enumerate(phonon_energies)
            if abs((ϵi-ϵf-ωph)*histogram_width)<0.5
                self_energy +=  π*abs(phonon_mat_elements[idx])^2*(1-fermi)*histogram_width/mesh^2
            end
            if abs((ϵi-ϵf+ωph)*histogram_width)<0.5
                self_energy +=  π*abs(phonon_mat_elements[idx])^2*(fermi)*histogram_width/mesh^2
            end
        end
    end
    return self_energy
end

"""
$(TYPEDSIGNATURES)

Find the phonon polarization for a 2d material- basically figure out which phonon branches correspond to logitudinal, transverse and z polarized modes
"""
function phonon_polarization(q::Vector{<:Float64}, eigenvector::Vector{<:ComplexF64}, lattice_vectors::Vector{<:Vector{<:Float64}})
    longitudinal_direction = unnormalize_kvector(lattice_vectors, q)
    lnorm = sqrt(sum(longitudinal_direction.^2))
    transverse_direction = [0 -1 0; 1 0 0; 0 0 1]*longitudinal_direction
    tnorm = lnorm
    !isequal(mod(length(eigenvector), 3), 0) && error("The eigenvector length must be a multiple of 3")
    natoms = div(length(eigenvector), 3)
    overlaps = [0,0,0]
    for i in 0:natoms-1
        eignorm = sqrt(sum((abs.(eigenvector[1+i*3:3+i*3])).^2))
        l = sum(eigenvector[1+i*3:3+i*3].*longitudinal_direction)/(lnorm*eignorm)
        t = sum(eigenvector[1+i*3:3+i*3].*transverse_direction)/(tnorm*eignorm)
        z = sum(eigenvector[1+i*3:3+i*3].*[0, 0, 1])/(eignorm)
        overlaps += abs.([l, t, z])
        println(abs(l)^2+abs(t)^2+abs(z)^2)
    end
    return overlaps
end