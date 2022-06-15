"""
$(TYPEDSIGNATURES)
"""
function export_hwannier(filebase::AbstractString, kmesh::Vector{<:Real}; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    if spin isa Val{'n'}
        cell_map = "$filebase.mlwfCellMap"
        cell_weights = "$filebase.mlwfCellWeights"
        H = "$filebase.mlwfH"
        band_file = "$filebase.txt"
        cell_map_file = "$filebase.map.txt"
    elseif spin isa Val{'u'}
        cell_map = "$filebase.mlwfCellMapUp"
        cell_weights = "$filebase.mlwfCellWeightsUp"
        H = "$filebase.mlwfHUp"
        band_file = "$(filebase)Up.txt"
        cell_map_file = "$(filebase)Up.map.txt"
    elseif spin isa Val{'d'}
        cell_map = "$filebase.mlwfCellMapDn"
        cell_weights = "$filebase.mlwfCellWeightsDn"
        H = "$filebase.mlwfHDn"
        band_file = "$(filebase)Dn.txt"
        cell_map_file = "$(filebase)Dn.map.txt"
    end

    cellMap = Int.(np.loadtxt(cell_map)[:,1:3])
    Wwannier = np.fromfile(cell_weights)
    nCells = first(size(cellMap))
    nBands = Int(sqrt(length(Wwannier) / nCells))
    Wwannier = reshape(permutedims(reshape(Wwannier, (nBands*nBands, nCells)), (2, 1)), (nCells, nBands, nBands))
    kfoldProd = prod(kmesh)
    kStride = [kmesh[2]*kmesh[3], kmesh[3], 1]
    Hreduced = np.fromfile(H)
    Hreduced = reshape(permutedims(reshape(Hreduced, (nBands*nBands, kfoldProd)), (2, 1)), (kfoldProd, nBands, nBands))
    iReduced = np.mod(cellMap, kmesh)*kStride
    Hwannier = Wwannier .* Hreduced[iReduced .+ 1, :, :]
    np.savetxt(band_file, np.reshape(Hwannier, (length(iReduced), nBands*nBands)))
    np.savetxt(cell_map_file, cellMap)
    return Hwannier, cellMap
end

export_hwannier(filebase; kwargs...) = export_hwannier(filebase, kpoint_folding("$filebase.out"); kwargs...)

"""
$(TYPEDSIGNATURES)
Exports momentum file in format filebase.momentum.txt

"""
function export_momentum(filebase::AbstractString, kmesh::Vector{<:Integer}; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    
    if spin isa Val{'n'}
        cell_map = filebase*".mlwfCellMap"
        cell_weights = filebase*".mlwfCellWeights"
        H = filebase*".mlwfH"
        P = filebase*".mlwfP"
        momentum_file = filebase*".momentum.txt"
    elseif spin isa Val{'u'}
        cell_map = filebase*".mlwfCellMapUp"
        cell_weights = filebase*".mlwfCellWeightsUp"
        H = filebase*".mlwfHUp"
        P = filebase*".mlwfPUp"
        momentum_file = filebase*"Up.momentum.txt"
    elseif spin isa Val{'d'}
        cell_map = filebase*".mlwfCellMapDn"
        cell_weights = filebase*".mlwfCellWeightsDn"
        H = filebase*".mlwfHDn"
        P = filebase*".mlwfPDn"    
        momentum_file = filebase*"Dn.momentum.txt"
    end
    cellMap = Int.(np.loadtxt(cell_map)[:,1:3])
    Wwannier = np.fromfile(cell_weights)
    nCells = first(size(cellMap))
    nBands = Int(sqrt(length(Wwannier) / nCells))
    Wwannier = reshape(permutedims(reshape(Wwannier, (nBands*nBands, nCells)), (2, 1)), (nCells, nBands, nBands))
    kfoldProd = prod(kmesh)
    kStride = [kmesh[2]*kmesh[3], kmesh[3], 1]
    iReduced = np.mod(cellMap, kmesh)*kStride
    Preduced = np.fromfile(P)
    Preduced = permutedims(reshape(np.fromfile(P), (3*nBands*nBands, kfoldProd)), (2, 1))
    Preduced = permutedims(reshape(Preduced, (kfoldProd, nBands, 3*nBands)), (1, 3, 2)) 
    Preduced = permutedims(reshape(Preduced, (kfoldProd, nBands, 3, nBands)), (1, 3, 4, 2 ))
    Pwannier = np.einsum("ijk, izjk -> izjk", Wwannier, Preduced[iReduced .+ 1, :,  :, :])
    np.savetxt(momentum_file, np.reshape(Pwannier, (length(iReduced), 3*nBands*nBands)))
    return Pwannier
end

export_momentum(filebase; kwargs...) = export_momentum(filebase, kpoint_folding("$filebase.out"); kwargs...)

"""
$(TYPEDSIGNATURES)

Exports heph, cellmapeph in files with pattern filebase.heph.txt and filebase.mapeph.txt
"""
function export_heph(filebase::AbstractString, qmesh::Vector{<:Integer}; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    
    if spin isa Val{'n'}
        cell_map = "$filebase.mlwfCellMap"
        cell_weights = "$filebase.mlwfCellWeights"
        cell_map_ph = "$filebase.mlwfCellMapPh"
        cell_map_ph_weights = "$filebase.mlwfCellWeightsPh"
        HPh = "$filebase.mlwfHePh"
        heph_file = "$filebase.heph.txt"
        celleph_file = "$filebase.mapeph.txt"
    elseif spin isa Val{'u'}
        cell_map = "$filebase.mlwfCellMapUp"
        cell_weights = "$filebase.mlwfCellWeightsUp"
        cell_map_ph = "$filebase.mlwfCellMapPhUp"
        cell_map_ph_weights = "$filebase.mlwfCellWeightsPhUp"
        HPh = "$filebase.mlwfHePhUp"        
        heph_file = "$(filebase)Up.heph.txt"
        celleph_file = "$(filebase)Up.mapeph.txt"
    elseif spin isa Val{'d'}
        cell_map = "$filebase.mlwfCellMapDn"
        cell_weights = "$filebase.mlwfCellWeightsDn"
        cell_map_ph = "$filebase.mlwfCellMapPhDn"
        cell_map_ph_weights = "$filebase.mlwfCellWeightsPhDn"
        HPh = "$filebase.mlwfHePhDn"   
        heph_file = "$(filebase)Dn.heph.txt"
        celleph_file = "$(filebase)Dn.mapeph.txt"
    end

    cellMap = Int.(np.loadtxt(cell_map)[:,1:3])
    Wwannier = np.fromfile(cell_weights)
    nCells = first(size(cellMap))
    nBands = Int(sqrt(length(Wwannier) / nCells))
    cellMapEph = Int.(np.loadtxt(cell_map_ph, usecols=[0,1,2]))
    nCellsEph = first(size(cellMapEph))
    prodPhononSup = prod(qmesh)
    phononSupStride = [qmesh[3]*qmesh[2], qmesh[3], 1]
    cellWeightsEph = np.fromfile(cell_map_ph_weights)
    nAtoms = div(length(cellWeightsEph), nCellsEph*nBands)
    nModes = 3*nAtoms
    cellWeightsEph = reshape(permutedims(reshape(cellWeightsEph, (nBands*nAtoms, nCellsEph)), (2, 1)), (nCellsEph, nAtoms, nBands))

    cellWeightsEph = np.repeat(reshape(cellWeightsEph, (nCellsEph,nAtoms,1,nBands)), 3, axis=2) #repeat atom weights for 3 directions
    
    cellWeightsEph = reshape(permutedims(cellWeightsEph, (1, 3, 2, 4)), (nCellsEph, nModes,nBands)) #coombine nAtoms x 3 into single dimension: nModes
        
    iReducedEph = np.mod(cellMapEph, qmesh)*phononSupStride

    HePhReduced = np.fromfile(HPh)

    HePhReduced = permutedims(reshape(HePhReduced, (nModes*nBands^2, prodPhononSup^2)), (2, 1))
    HePhReduced = permutedims(reshape(HePhReduced, (prodPhononSup, prodPhononSup, nModes*nBands^2)), (2, 1, 3))
    HePhReduced = permutedims(reshape(HePhReduced, (prodPhononSup, prodPhononSup, nBands^2, nModes)), (1, 2, 4, 3))
    HePhReduced = reshape(HePhReduced, (prodPhononSup, prodPhononSup, nModes, nBands, nBands))

    HePhWannier = np.einsum("ikl, jkm -> ijklm", cellWeightsEph, cellWeightsEph)
    HePhWannier = HePhWannier .* HePhReduced[iReducedEph  .+ 1, :, :, :, :][:, iReducedEph .+ 1, :, :, :]
    np.savetxt(heph_file, np.reshape(HePhWannier, (length(iReducedEph)^2, -1)))
    np.savetxt(celleph_file, cellMapEph)

    return HePhWannier, cellMapEph
end

export_heph(filebase; kwargs...) = export_heph(filebase, phonon_supercell("$filebase.out"); kwargs...)
