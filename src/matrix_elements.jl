function eph_matrix_elements(Heph::Array{<:Real, 5}, celleph_map::Matrix{<:Real}, U1::Matrix{<:ComplexF64}, U2::Matrix{<:ComplexF64}, 
    omegaph::Vector{<:Float64}, Uph::Matrix{<:ComplexF64}, k1::Vector{<:Real}, k2::Vector{<:Real})
    phase1 = exp.((2im*π)*(celleph_map*k1))
    phase2 = exp.((2im*π)*(celleph_map*k2))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaph, 1e-6))

    g = np.einsum("bd, ycb-> ycd", U2, 
    np.einsum("ac,yab -> ycb", conj(U1),
    np.einsum("xy, xab-> yab", Uph, 
    np.einsum("R, Rxab -> xab", phase2, 
    np.einsum("r, rRxab -> Rxab", conj(phase1), 
    Heph))))).*normFac 
    return g/eV
end

function eph_matrix_elements(Heph::Array{<:Real, 5}, celleph_map::Matrix{<:Real}, U1s::Array{<:ComplexF64, 3}, U2s::Array{<:ComplexF64, 3}, 
    omegaphs::Array{<:Float64, 3}, Uphs::Array{<:ComplexF64, 4}, k1s::AbstractMatrix{<:Real}, k2s::AbstractMatrix{<:Real})
    
    phase1 = exp.((2im*π)*(np.tensordot(celleph_map, k1s, axes=1)))
    phase2 = exp.((2im*π)*(np.tensordot(celleph_map, k2s, axes=1)))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaphs,1e-6))

    g = np.einsum("qbd, kqycb-> kqycd", U2s, 
    np.einsum("kac, kqyab -> kqycb", conj(U1s),
    np.einsum("kqxy, kqxab-> kqyab", Uphs, 
    np.einsum("Rq, kRxab -> kqxab", phase2, 
    np.einsum("rk, rRxab -> kRxab", conj(phase1), 
    Heph))))) .* normFac 

    return g/eV
end

function eph_matrix_elements(Heph::Array{<:Real, 5}, celleph_map::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, 
    cellph_map::Array{<:Real, 2}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    k1::Vector{<:Real}, k2::Vector{<:Real})

    _, U2 = wannier_bands(Hwannier, cell_map, k2) 
    _, U1 = wannier_bands(Hwannier, cell_map, k1) 

    omegaphsquared, Uph = diagonalize_phonon(force_matrix, cellph_map, k1-k2)
    omegaph = sqrt.(abs.(omegaphsquared))
    return eph_matrix_elements(Heph, celleph_map, U1, U2, omegaph, Uph, k1, k2)
end

function eph_matrix_elements(Heph::Array{<:Real, 5}, celleph_map::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, 
    cellph_map::Array{<:Real, 2}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    k1s::AbstractMatrix{<:Real}, k2s::AbstractMatrix{<:Real})

    _, U2s = wannier_bands(Hwannier, cell_map, k2s) 
    _, U1s = wannier_bands(Hwannier, cell_map, k1s) 

    omegaphsquareds, Uphs = diagonalize_phonon(force_matrix, cellph_map, k1s, k2s)
    omegaphs = sqrt.(abs.(omegaphsquareds))
    return eph_matrix_elements(Heph, celleph_map, U1s, U2s, omegaphs, Uphs, k1s, k2s)
end

"""
$(TYPEDSIGNATURES)
"""
function momentum_matrix_elements(Pwannier::Array{Float64, 4}, cell_map::Array{Float64, 2}, k::Vector{<:Real})
    phase = np.exp(2im*π*cell_map*k); 
    Pk = np.tensordot(phase, Pwannier, axes=1); 
    return Pk*ħ/bohrtoangstrom
end

function momentum_from_bloch(lattice_vectors::Vector{<:Vector{<:Real}}, Hwannier::Array{Float64, 3}, cellmap::Matrix{Float64}, 
    k::Vector{<:Real}; qdiff::Vector{<:Real}=[1, 0, 0])
    
    qabs = norm(unnormalize_kvector(lattice_vectors, qdiff))
    prefactor = mₑ/ħ
    Energieskplusq1, Uskplusq1 = wannier_bands(Hwannier, cellmap, k + qdiff/2)
    Energieskplusq2, Uskplusq2 = wannier_bands(Hwannier, cellmap, k - qdiff/2)

    numbands = size(Hwannier)[2]
    EnergiesDiff = np.repeat(np.reshape(Energieskplusq1, (numbands, 1)), numbands, axis = 1) -
                np.repeat(np.reshape(Energieskplusq2, (1, numbands)), numbands, axis = 0)
    p = prefactor/qabs * np.einsum("nm, nm -> nm", EnergiesDiff, np.einsum("in, im -> nm", np.conj(Uskplusq2), Uskplusq1))
    return p
end

function momentum_matrix_elements(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, k::Vector{<:Real})
    phase = np.exp(2im*π*cell_map*k); 
    Pk = np.tensordot(phase, Pwannier, axes=1); 
    _, Us = wannier_bands(Hwannier, cell_map, k) 
    Pk = np.einsum("ba, pbc, cd-> pad", conj(Us), Pk, Us)  
    return Pk*ħ/bohrtoangstrom
end

function momentum_matrix_elements(Us::Array{<:ComplexF64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, kpoints::AbstractArray{<:Real, 2})
    phase = np.exp(2im*π*cell_map*kpoints); 
    Pk = np.einsum("al, aijk -> lijk", phase, Pwannier)
    Pk = np.einsum("kba, kpbc, kcd-> kpad", conj(Us), Pk, Us)  
    return Pk*ħ/bohrtoangstrom
end

function momentum_matrix_elements(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, kpoints::AbstractArray{<:Real, 2})
    _, Us = wannier_bands(Hwannier, cell_map, kpoints) 
    return momentum_matrix_elements(Us, cell_map, Pwannier, kpoints)
end

function pwannier(momentum_file::AbstractString, cell_map_file::AbstractString) 
    cell_map_numlines = countlines(cell_map_file);
    numbands = Int(sqrt(np.size(np.loadtxt(momentum_file))//(cell_map_numlines*3)));
    Pwannier = np.reshape(np.loadtxt(momentum_file), (cell_map_numlines, 3, numbands, numbands));
    return Pwannier
end

pwannier(filebase::AbstractString) = pwannier("$filebase.momentum.txt", "$filebase.map.txt")

"""
$(TYPEDSIGNATURES)
"""
function hephwannier(filebase::AbstractString) 
    heph_file = "$filebase.heph.txt"
    celleph_file = "$filebase.mapeph.txt"
    wannier_file = "$filebase.txt"
    cell_map_file = "$filebase.map.txt"
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = np.loadtxt(wannier_file)
    numbands =  length(Hwannier) == cell_map_numlines ? 1 : Int(sqrt(size(Hwannier)[2]))
    celleph_map_numlines = countlines(celleph_file)
    Heph = np.reshape(np.loadtxt(heph_file), (celleph_map_numlines, celleph_map_numlines, -1, numbands, numbands))
    return Heph
end

