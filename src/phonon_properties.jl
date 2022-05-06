"""
$(TYPEDSIGNATURES)
Plots the phonon band dispersion at the kpoints supplied
"""
function plot_phonons(cell_map::AbstractString, phononOmegaSq::AbstractString, kpoints::AbstractString="bandstruct.kpoints"; kwargs...)
    cellMapPh = np.loadtxt(cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    kpointsIn = np.loadtxt(kpoints, skiprows=2, usecols=(1,2,3))
    forceMatrixTilde = np.tensordot(np.exp((2im*np.pi)*np.dot(kpointsIn,transpose(cellMapPh))), forceMatrixPh, axes=1)
    omegaSq, _ = np.linalg.eigh(forceMatrixTilde)
    plot(title="Phonon Dispersion", titlefontsize=20, ytickfontsize=15,  yguidefontsize=30, sqrt.(abs.(omegaSq))/eV, ylabel= "Energy (eV)", linewidth=2, color="orange", legend=false, size=(800, 1000), xticks=[]; kwargs...)
end

"""
$(TYPEDSIGNATURES)
Give phonon dispersion at individual kpoints
"""
function phonon_dispersion(phonon_cell_map::AbstractString, phononOmegaSq::AbstractString, qnorm::Vector{<:Real}) 
    cellMapPh = np.loadtxt(phonon_cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    #--- Fourier transform from real to k space:
    forceMatrixTildeq = np.tensordot(np.exp((2im*np.pi)*np.dot(qnorm,transpose(cellMapPh))), forceMatrixPh, axes=1)
    omegaSq, _ = np.linalg.eigh(forceMatrixTildeq)
    return sqrt.(abs.(omegaSq))/eV
end

"""
$(TYPEDSIGNATURES)

Optional keyword argument return_negative may be supplied to return imaginary frequencies as negative frequencies

"""
function phonon_dispersion(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorm::Vector{<:Real}; return_negative::Bool=false)
    forceMatrixTildeq = np.tensordot(np.exp(2im*π*np.dot(qnorm, transpose(phonon_cell_map))), force_matrix, axes=1)
    omegaSq, _ = np.linalg.eigh(forceMatrixTildeq)
    freq = return_negative ? sign.(omegaSq).*sqrt.(abs.(omegaSq))/eV : sqrt.(abs.(omegaSq))/eV
    return freq
end

"""
$(TYPEDSIGNATURES)
Returns phonon dispersion along a supplied q-path
"""
function phonon_dispersionpath(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorms::Vector{<:Vector{<:Real}}; return_negative=false)
    nmodes = length(phonon_dispersion(force_matrix, phonon_cell_map, [0, 0, 0]))
    nks = length(qnorms)
    ϵalongpath = zeros(nks, nmodes )
    for (index, qnorm) in enumerate(qnorms)
        ϵalongpath[index, :] = phonon_dispersion(force_matrix, phonon_cell_map, qnorm, return_negative=return_negative)
    end
    return ϵalongpath
end

"""
$(TYPEDSIGNATURES)

Plots the phonon dispersion along a (possibly interpolated) Brillouin zone path as read from a jdftx convention bandstruct.kpoints file.
"""
function phonon_dispersionpath(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}; kpointsfile::AbstractString="bandstruct.kpoints", 
    interpolate::Integer=1, plotbands::Bool=false, return_negative=false, kwargs...)
    
    nmodes = length(phonon_dispersion(force_matrix, phonon_cell_map, [0, 0, 0]))
    qnorms = bandstructkpoints2q(kpointsfile=kpointsfile, interpolate=interpolate)
    nks = length(qnorms)
    ϵalongpath = zeros(nks, nmodes )
    for (index, qnorm) in enumerate(qnorms)
        ϵalongpath[index, :] = phonon_dispersion(force_matrix, phonon_cell_map, qnorm, return_negative=return_negative)
    end
    plotbands && plot(ϵalongpath; kwargs...)
    return ϵalongpath
end

"""
$(TYPEDSIGNATURES)
"""
function phonon_dispersionpath(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorms::Tuple{Vector{<:Vector{<:Real}}, Integer})
    diffs = diff(qnorms[1])/qnorms[2]
    allqnorms = Vector{Vector{Float64}}()
    println(qnorms[2])
    for (diff, k) in zip(diffs, qnorms[1][1:end-1])
        for interp in 0:qnorms[2]-1
            push!(allqnorms, k+interp*diff)
        end
    end
    return phonon_dispersionpath(force_matrix, phonon_cell_map, allqnorms)
end

"""
$(TYPEDSIGNATURES)
"""
function phonon_dispersionmodes(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorm::Vector{<:Real})
    phase = np.exp((2im*np.pi)*np.tensordot(qnorm, transpose(phonon_cell_map), axes=1))
    ### Note that we must permute the indices of the force matrix by jdftx convention
    ### As in http://jdftx.org/EphMatrixElements.html
    omegaSq, U = np.linalg.eigh(np.tensordot(phase, permutedims(force_matrix, (1, 3, 2)), axes=1))
    return sqrt.(abs.(omegaSq))/eV, U
end

"""
$(TYPEDSIGNATURES)

"""
function phonon_force_matrix(phonon_cell_map::AbstractString, phononOmegaSq::AbstractString)
    cellMapPh = np.loadtxt(phonon_cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    return forceMatrixPh, cellMapPh
end

"""
$(TYPEDSIGNATURES)
"""
function phonon_force_matrix(filebase::AbstractString)
    phonon_cell_map = "$(filebase).phononCellMap"
    phononOmegaSq = "$(filebase).phononOmegaSq"
    cellMapPh = np.loadtxt(phonon_cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))
    println("Number of phonon modes is: ", nModesPh, "\nIf this is incorrect, something went wrong somewhere at some point.")
    return forceMatrixPh, cellMapPh
end


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