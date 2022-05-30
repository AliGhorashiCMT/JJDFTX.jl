"""
$(TYPEDSIGNATURES)
"""
function write_map_write_h(filebase::AbstractString, kmesh::Vector{<:Real}; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    if spin isa Val{'u'}
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
    elseif spin isa Val{'n'}
        cell_map = "$filebase.mlwfCellMap"
        cell_weights = "$filebase.mlwfCellWeights"
        H = "$filebase.mlwfH"
        band_file = "$filebase.txt"
        cell_map_file = "$filebase.map.txt"
    end
    py"""   
    def write_map_write_h_py(cell_map, cell_weights, H, kmesh, band_file, cell_map_file):
        import numpy as np
        cellMapUp = np.loadtxt(cell_map)[:,0:3].astype(np.int)
        WwannierUp = np.fromfile(cell_weights)
        nCellsUp = cellMapUp.shape[0]
        nBandsUp = int(np.sqrt(WwannierUp.shape[0] / nCellsUp))
        WwannierUp = WwannierUp.reshape((nCellsUp,nBandsUp,nBandsUp)).swapaxes(1,2)
        kfold=np.array([kmesh[0], kmesh[1], kmesh[2]])
        kfoldProd = np.prod(kfold)
        kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
        HreducedUp = np.fromfile(H).reshape((kfoldProd,nBandsUp,nBandsUp)).swapaxes(1,2)
        iReducedUp = np.dot(np.mod(cellMapUp, kfold[None,:]), kStride)
        HwannierUp = WwannierUp * HreducedUp[iReducedUp]
        np.savetxt(band_file, HwannierUp.reshape(len(iReducedUp), nBandsUp*nBandsUp ))
        np.savetxt(cell_map_file, cellMapUp)
    """
    py"write_map_write_h_py"(cell_map, cell_weights, H, kmesh, band_file, cell_map_file)
end

"""
$(TYPEDSIGNATURES)
"""
function export_hwannier(filebase::AbstractString, kmesh::Vector{<:Real}; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    if spin isa Val{'u'}
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
    elseif spin isa Val{'n'}
        cell_map = "$filebase.mlwfCellMap"
        cell_weights = "$filebase.mlwfCellWeights"
        H = "$filebase.mlwfH"
        band_file = "$filebase.txt"
        cell_map_file = "$filebase.map.txt"
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
end


"""
$(TYPEDSIGNATURES)

"""
function write_momentum(filebase::AbstractString, kmesh::Vector{<:Integer}, momentum_file::AbstractString; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    
    if spin isa Val{'u'}
        cell_map = filebase*".mlwfCellMap"
        cell_weights = filebase*".mlwfCellWeights"
        H = filebase*".mlwfH"
        P = filebase*".mlwfP"
    elseif spin isa Val{'u'}
        cell_map = filebase*".mlwfCellMapUp"
        cell_weights = filebase*".mlwfCellWeightsUp"
        H = filebase*".mlwfHUp"
        P = filebase*".mlwfPUp"
    elseif spin isa Val{'d'}
        cell_map = filebase*".mlwfCellMapDn"
        cell_weights = filebase*".mlwfCellWeightsDn"
        H = filebase*".mlwfHDn"
        P = filebase*".mlwfPDn"    
    end
    py"""
    def write_map_write_p(cell_map, cell_weights, H, P, kmesh, momentum_file):
        import numpy as np
        cellMap = np.loadtxt(cell_map)[:,0:3].astype(np.int)
        Wwannier = np.fromfile(cell_weights)
        nCells = cellMap.shape[0]
        nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
        Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)
        kfold=np.array([kmesh[0], kmesh[1], kmesh[2]])
        kfoldProd = np.prod(kfold)
        kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
        Hreduced = np.fromfile(H).reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)
        Preduced = np.fromfile(P).reshape((kfoldProd,3,nBands,nBands)).swapaxes(2,3)
        iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
        Hwannier = Wwannier * Hreduced[iReduced]
        Pwannier = Wwannier[:,None] * Preduced[iReduced]
        np.savetxt(momentum_file, Pwannier.reshape(len(iReduced), 3*nBands*nBands ))
    """
    py"write_map_write_p"(cell_map, cell_weights, H, P, kmesh, momentum_file)
end

"""
$(TYPEDSIGNATURES)
"""
function write_eph_matrix_elements(filebase::AbstractString, nModes::Integer, qmesh::Vector{<:Integer}; 
    spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    
    if spin isa Val{'u'}
        cell_map = "$filebase.mlwfCellMap"
        cell_weights = "$filebase.mlwfCellWeights"
        cell_map_ph = "$filebase.mlwfCellMapPh"
        cell_map_ph_weights = "$filebase.mlwfCellWeightsPh"
        HPh = "$filebase.mlwfHePh"
    elseif spin isa Val{'u'}
        cell_map = "$filebase.mlwfCellMapUp"
        cell_weights = "$filebase.mlwfCellWeightsUp"
        cell_map_ph = "$filebase.mlwfCellMapPhUp"
        cell_map_ph_weights = "$filebase.mlwfCellWeightsPhUp"
        HPh = "$filebase.mlwfHePhUp"        
    elseif spin isa Val{'d'}
        cell_map = "$filebase.mlwfCellMapDn"
        cell_weights = "$filebase.mlwfCellWeightsDn"
        cell_map_ph = "$filebase.mlwfCellMapPhDn"
        cell_map_ph_weights = "$filebase.mlwfCellWeightsPhDn"
        HPh = "$filebase.mlwfHePhDn"    
    end
    py"""
    def write_eph(cell_map, cell_weights, cell_map_ph, cell_map_ph_weights, HPh, nModes, qmesh):
        import numpy as np
        cellMap = np.loadtxt(cell_map)[:,0:3].astype(np.int)
        Wwannier = np.fromfile(cell_weights)
        nCells = cellMap.shape[0]
        nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
        cellMapEph = np.loadtxt(cell_map_ph, usecols=[0,1,2]).astype(int)
        nCellsEph = cellMapEph.shape[0]
        prodPhononSup = np.prod(qmesh)
        phononSupStride = np.array([qmesh[1]*qmesh[2], qmesh[2], 1])
        nAtoms = nModes // 3 # Rounding down to an integer, since we know the number of modes is a multiple of three
        cellWeightsEph = np.fromfile(cell_map_ph_weights).reshape((nCellsEph,nBands,nAtoms)).swapaxes(1,2)
        cellWeightsEph = np.repeat(cellWeightsEph.reshape((nCellsEph,nAtoms,1,nBands)), 3, axis=2) #repeat atom weights for 3 directions
        cellWeightsEph = cellWeightsEph.reshape((nCellsEph,nModes,nBands)) #coombine nAtoms x 3 into single dimension: nModes
        iReducedEph = np.dot(np.mod(cellMapEph, qmesh[None,:]), phononSupStride)
        HePhReduced = np.fromfile(HPh).reshape((prodPhononSup,prodPhononSup,nModes,nBands,nBands)).swapaxes(3,4)
        HePhWannier = cellWeightsEph[:,None,:,:,None] * cellWeightsEph[None,:,:,None,:] * HePhReduced[iReducedEph][:,iReducedEph]
        return HePhWannier, cellMapEph
    """
    py"write_eph"(cell_map, cell_weights, cell_map_ph, cell_map_ph_weights, HPh, nModes, qmesh)
end
