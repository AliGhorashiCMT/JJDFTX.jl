#function phmatrixelements(k1::Array{T, 1}, k2::Array{R, 1}) where {T<:Number, R<:Number}
function phmatrixelements(k1, k2) 
    py"""
    import numpy as np
    cellMap = np.loadtxt("wannierDefect.mlwfCellMapUp")[:,0:3].astype(int)
    Wwannier = np.fromfile("wannierDefect.mlwfCellWeightsUp")
    nCells = cellMap.shape[0]
    nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
    ##print(nBands)
    Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)
    #--- Get k-point folding from totalE.out:
    kfold = np.array([6, 6, 1])
    kfoldProd = np.prod(kfold)
    kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
    #--- Read reduced Wannier Hamiltonian, momenta and expand them:
    Hreduced = np.fromfile("wannierDefect.mlwfHUp").reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)
    iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
    Hwannier = Wwannier * Hreduced[iReduced]
    #Read phonon dispersion relation:
    cellMapPh = np.loadtxt('BN33BC.phononCellMap', usecols=[0,1,2]).astype(int)
    nCellsPh = cellMapPh.shape[0]
    omegaSqR = np.fromfile('BN33BC.phononOmegaSq') #just a list of numbers
    nModes = int(np.sqrt(omegaSqR.shape[0] // nCellsPh))
    omegaSqR = omegaSqR.reshape((nCellsPh, nModes, nModes)).swapaxes(1,2)
    #Read e-ph matrix elements
    cellMapEph = np.loadtxt('wannierDefect.mlwfCellMapPhUp', usecols=[0,1,2]).astype(int)
    nCellsEph = cellMapEph.shape[0]
    #--- Get phonon supercell from phonon.out:
    phononSup = np.array([1, 1, 1])
    prodPhononSup = np.prod(phononSup)
    phononSupStride = np.array([phononSup[1]*phononSup[2], phononSup[2], 1])
    #--- Read e-ph cell weights:
    nAtoms = nModes // 3
    ##print(nAtoms)
    ##print(nCellsEph)
    cellWeightsEph = np.fromfile("wannierDefect.mlwfCellWeightsPhUp").reshape((nCellsEph,nBands,nAtoms)).swapaxes(1,2)
    cellWeightsEph = np.repeat(cellWeightsEph.reshape((nCellsEph,nAtoms,1,nBands)), 3, axis=2) #repeat atom weights for 3 directions
    cellWeightsEph = cellWeightsEph.reshape((nCellsEph,nModes,nBands)) #coombine nAtoms x 3 into single dimension: nModes
    #--- Read, reshape and expand e-ph matrix elements:
    iReducedEph = np.dot(np.mod(cellMapEph, phononSup[None,:]), phononSupStride)
    HePhReduced = np.fromfile('wannierDefect.mlwfHePhUp').reshape((prodPhononSup,prodPhononSup,nModes,nBands,nBands)).swapaxes(3,4)
    HePhWannier = cellWeightsEph[:,None,:,:,None] * cellWeightsEph[None,:,:,None,:] * HePhReduced[iReducedEph][:,iReducedEph]
    '''
    print(np.shape(HePhWannier))
    '''
    #Constants / calculation parameters:
    eV = 1/27.2114 #in Hartrees
    #Calculate energies, eigenvectors and velocities for given k
    def calcE(k):
        #Fourier transform to k:
        phase = np.exp((2j*np.pi)*np.dot(k,cellMap.T))
        H = np.tensordot(phase, Hwannier, axes=1)
        #Diagonalize and switch to eigen-basis:
        E,U = np.linalg.eigh(H) #Diagonalize
        return E, U
    #Calculate phonon energies and eigenvectors for given q
    def calcPh(q):
        phase = np.exp((2j*np.pi)*np.tensordot(q,cellMapPh.T, axes=1))
        omegaSq, U = np.linalg.eigh(np.tensordot(phase, omegaSqR, axes=1))
        omegaPh = np.sqrt(np.maximum(omegaSq, 0.))
        return omegaPh, U
    #Calculate e-ph matrix elements, along with ph and e energies, and e velocities
    def calcEph(k1, k2):
        #Electrons:
        #E1, U1= calcE(k1)
        #E2, U2 = calcE(k2)
        #print(U1)
        #print(U2)
        #Phonons for all pairs pf k1 - k2:
        #omegaPh, Uph = calcPh(k1[:,None,:] - k2[None,:,:])
        omegaPh, Uph = calcPh(k1-k2)
        '''
        print(omegaPh)
        print(np.shape(Uph))
        '''
        #E-ph matrix elements for all pairs of k1 - k2:
        phase1 = np.exp((2j*np.pi)*np.dot(k1,cellMapEph.T))
        phase2 = np.exp((2j*np.pi)*np.dot(k2,cellMapEph.T))

        '''
        print(phase1)
        print(phase2)
        '''
        normFac = np.sqrt(0.5/np.maximum(omegaPh,1e-6))
        '''
        print(normFac)
        '''
        '''
        g1 = np.einsum('Kbd,kKycb->kKycd', U2, #Rotate to electron 2 eigenbasis
            np.einsum('kac,kKyab->kKycb', U1.conj(), #Rotate to electron 1 eigenbasis
            np.einsum('kKxy,kKxab->kKyab', Uph, #Rotate to phonon eigenbasis
            np.einsum('KR,kRxab->kKxab', phase2, #Fourier transform from r2 -> k2
            np.einsum('kr,rRxab->kRxab', phase1.conj(), #Fourier transform from r1 -> k1
            HePhWannier)))))*normFac[...,None,None] #Phonon amplitude factor
        '''


        '''
        The following g is from only diagonalizing in the phonon basis and performing a fourier transform. This is appropriate
        if there is only one wannier band in question
        '''
        
        '''
        g2 = np.einsum('kKxy,kKxab->kKyab', Uph, np.einsum('KR,kRxab->kKxab', phase2, np.einsum('kr,rRxab->kRxab', phase1.conj(),  HePhWannier)))*normFac[...,None,None] #Phonon amplitude factor
        '''

        '''
            We alter the einsum such that the convention is that the kvectors given are three element arrays instead of arrays of 
            three element arrays.
        '''
        g2 = np.einsum('xy, xab-> yab', Uph, np.einsum('R,Rxab->xab', phase2, np.einsum('r,rRxab->Rxab', phase1.conj(),  HePhWannier)))

        #return g2
        #return Uph
        return g2.flatten()/eV*normFac
        #return HePhWannier
        #return g1, g2, omegaPh
    """
    py"calcEph"(k1, k2) 

end

#= 
Next, we will examine the momentum matrix elements. 
Note that the matrix elements are initially in the Wannier basis and must be transformed to the Bloch basis!
=#
function eph_matrix_elements(HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, k1::Array{<:Real, 1}, k2::Array{<:Real, 1})
    #Phonons for all pairs pf k1 - k2:
    #omegaPh, Uph = calcPh(k1[:,None,:] - k2[None,:,:])
    omegaPh, Uph = phonon_dispersionmodes(force_matrix, phonon_cell_map, k1-k2)
    ##Note that the phonon energies given by phonon dispersionmodes are in eV, so they must be converted 
    omegaPh *= eV
    #phase1 = np.exp((2im*π )*(cellMapEph*k1))
    phase1 = exp.((2im*π )*(cellMapEph*k1))
    #phase2 = np.exp((2im*π)*(cellMapEph*k2))
    #phase1 = np.exp((2im*np.pi)*np.dot(k1, transpose(cellMapEph)))
    #phase2 = np.exp((2im*np.pi)*np.dot(k2, transpose(cellMapEph)))
    phase2 = exp.((2im*π)*(cellMapEph*k2))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaPh,1e-6))
    #normFac = np.sqrt(0.5/max.(omegaPh, Ref(1e-6)))
    #print(typeof(normFac))
    #=g=  np.einsum("kKxy,kKxab->kKyab", Uph, #Rotate to phonon eigenbasis
        np.einsum("KR,kRxab->kKxab", phase2, #Fourier transform from r2 -> k2
        np.einsum("kr,rRxab->kRxab", phase1.conj(), #Fourier transform from r1 -> k1
        HePhWannier))) * normFac #Phonon amplitude factor
    =#
    g= vec(np.einsum("xy, xab-> yab", Uph, #Rotate to phonon eigenbasis
        np.einsum("R,Rxab->xab", phase2, #Fourier transform from r2 -> k2
        np.einsum("r,rRxab->Rxab", conj(phase1), #Fourier transform from r1 -> k1
        HePhWannier)))).*normFac  #Phonon amplitude factor
    return g/eV
end

function eph_matrix_elements(HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, wannier_file::String, cell_map_file::String, k1::Array{<:Real, 1}, k2::Array{<:Real, 1}, nbands::Int)
    omegaPh, Uph = phonon_dispersionmodes(force_matrix, phonon_cell_map, k1-k2)
    ##Note that the phonon energies given by phonon dispersionmodes are in eV, so they must be converted 
    omegaPh *= eV
    phase1 = exp.((2im*π )*(cellMapEph*k1))
    phase2 = exp.((2im*π)*(cellMapEph*k2))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaPh,1e-6))
    U2 = wannier_vectors(wannier_file::String, cell_map_file::String, k2, nbands) 
    U1 = wannier_vectors(wannier_file::String, cell_map_file::String, k2, nbands) 
    g = np.einsum("bd, ycb-> ycd", U2, #Rotate to electron 2 eigenbasis
    np.einsum("ac,yab -> ycb", conj(U1), #Rotate to electron 1 eigenbasis
    np.einsum("xy, xab-> yab", Uph, #Rotate to phonon eigenbasis
    np.einsum("R, Rxab -> xab", phase2, #Fourier transform from r2 -> k2
    np.einsum("r, rRxab -> Rxab", conj(phase1), #Fourier transform from r1 -> k1
    HePhWannier))))).*normFac #Phonon amplitude factor
    #=    
    g= vec(np.einsum("xy, xab-> yab", Uph, #Rotate to phonon eigenbasis
        np.einsum("R,Rxab->xab", phase2, #Fourier transform from r2 -> k2
        np.einsum("r,rRxab->Rxab", conj(phase1), #Fourier transform from r1 -> k1
        HePhWannier)))).*normFac  #Phonon amplitude factor
    =#
    return g/eV
end

function momentum_matrix_elements(Pwannier::Array{Float64, 4}, cell_map::Array{Float64, 2}, k::Array{<:Real, 1})
    phase = np.exp(2im*π*cell_map*k); 
    Pk = np.tensordot(phase, Pwannier, axes=1); 
    #= 
    JDFTX output is in Atomic units. Therefore, the units of the momentum matrix elements are in ħ/a₀
    a₀ is the Bohr radius, which is approximately 0.529 Angstrom. ħ is 6.6*10-16 eV*seconds. Since all 
    physical calculations perfomed in jdftx_to_plot assume lengths to be given in angstroms and energies 
    to be given in eV, we multiply Pk by ħ/bohrtoangstrom
    =#
    return Pk*ħ/bohrtoangstrom
end

function momentum_matrix_elements(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, k::Array{<:Real, 1})
    phase = np.exp(2im*π*cell_map*k); 
    Pk = np.tensordot(phase, Pwannier, axes=1); 
    Vk = wannier_vectors(Hwannier, cell_map, k) 
    Pk = np.einsum("ba, pbc, cd-> pad",   #Sum using Einstein notation for
    conj(Vk), Pk, Vk)  
    #= 
    JDFTX output is in Atomic units. Therefore, the units of the momentum matrix elements are in ħ/a₀
    a₀ is the Bohr radius, which is approximately 0.529 Angstrom. ħ is 6.6*10-16 eV*seconds. Since all 
    physical calculations perfomed in jdftx_to_plot assume lengths to be given in angstroms and energies 
    to be given in eV, we multiply Pk by ħ/bohrtoangstrom
    =#
    return Pk*ħ/bohrtoangstrom
end

function pwannier(pwannier_file::String, cell_map_file::String, nbands::Int64) 
    cell_map = np.loadtxt(cell_map_file)
    cell_map_numlines = countlines(cell_map_file)
    Pwannier = np.reshape(np.loadtxt(pwannier_file), (cell_map_numlines, 3, nbands, nbands))
    return Pwannier
end
