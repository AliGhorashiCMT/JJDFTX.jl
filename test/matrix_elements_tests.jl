@testset "Checking Electron Phonon Matrix Elements with Shankar's Implementation in Python" begin
    py"""
    import numpy as np
    dir = "../data/electron_phonon_reference/"
    cellMap = np.loadtxt(dir+"wannier.mlwfCellMap")[:,0:3].astype(int)
    Wwannier = np.fromfile(dir+"wannier.mlwfCellWeights")
    nCells = cellMap.shape[0]
    nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
    Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)
    for line in open(dir+'totalE.out'):
        if line.startswith('kpoint-folding'):
            kfold = np.array([int(tok) for tok in line.split()[1:4]])
    kfoldProd = np.prod(kfold)
    kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
    Hreduced = np.fromfile(dir+"wannier.mlwfH").reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)
    Preduced = np.fromfile(dir+"wannier.mlwfP").reshape((kfoldProd,3,nBands,nBands)).swapaxes(2,3)
    iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
    Hwannier = Wwannier * Hreduced[iReduced]
    Pwannier = Wwannier[:,None] * Preduced[iReduced]
    cellMapPh = np.loadtxt(dir+'totalE.phononCellMap', usecols=[0,1,2]).astype(int)
    nCellsPh = cellMapPh.shape[0]
    omegaSqR = np.fromfile(dir+'totalE.phononOmegaSq') #just a list of numbers
    nModes = int(np.sqrt(omegaSqR.shape[0] // nCellsPh))
    omegaSqR = omegaSqR.reshape((nCellsPh, nModes, nModes)).swapaxes(1,2)
    cellMapEph = np.loadtxt(dir+'wannier.mlwfCellMapPh', usecols=[0,1,2]).astype(int)
    nCellsEph = cellMapEph.shape[0]
    phononSup = np.array([2, 2, 2])
    prodPhononSup = np.prod(phononSup)
    phononSupStride = np.array([phononSup[1]*phononSup[2], phononSup[2], 1])
    nAtoms = nModes // 3
    cellWeightsEph = np.fromfile(dir+"wannier.mlwfCellWeightsPh").reshape((nCellsEph,nBands,nAtoms)).swapaxes(1,2)
    cellWeightsEph = np.repeat(cellWeightsEph.reshape((nCellsEph,nAtoms,1,nBands)), 3, axis=2) #repeat atom weights for 3 directions
    cellWeightsEph = cellWeightsEph.reshape((nCellsEph,nModes,nBands)) #coombine nAtoms x 3 into single dimension: nModes
    iReducedEph = np.dot(np.mod(cellMapEph, phononSup[None,:]), phononSupStride)
    HePhReduced = np.fromfile(dir+'wannier.mlwfHePh').reshape((prodPhononSup,prodPhononSup,nModes,nBands,nBands)).swapaxes(3,4)
    HePhWannier = cellWeightsEph[:,None,:,:,None] * cellWeightsEph[None,:,:,None,:] * HePhReduced[iReducedEph][:,iReducedEph]
    def calcE(k):
        phase = np.exp((2j*np.pi)*np.dot(k,cellMap.T))
        H = np.tensordot(phase, Hwannier, axes=1)
        P = np.tensordot(phase, Pwannier,  axes=1)
        E,U = np.linalg.eigh(H) #Diagonalize
        v = np.imag(np.einsum('kba,kibc,kca->kai', U.conj(), P, U)) #diagonal only
        return E, U, v
    def calcPh(q):
        phase = np.exp((2j*np.pi)*np.tensordot(q,cellMapPh.T, axes=1))
        omegaSq, U = np.linalg.eigh(np.tensordot(phase, omegaSqR, axes=1))
        omegaPh = np.sqrt(np.maximum(omegaSq, 0.))
        return omegaPh, U
    def calcEph(k1, k2):
        E1, U1, v1 = calcE(k1)
        E2, U2, v2 = calcE(k2)
        omegaPh, Uph = calcPh(k1[:,None,:] - k2[None,:,:])
        phase1 = np.exp((2j*np.pi)*np.dot(k1,cellMapEph.T))
        phase2 = np.exp((2j*np.pi)*np.dot(k2,cellMapEph.T))
        normFac = np.sqrt(0.5/np.maximum(omegaPh,1e-6))
        g = np.einsum('Kbd,kKycb->kKycd', U2, #Rotate to electron 2 eigenbasis
            np.einsum('kac,kKyab->kKycb', U1.conj(), #Rotate to electron 1 eigenbasis
            np.einsum('kKxy,kKxab->kKyab', Uph, #Rotate to phonon eigenbasis
            np.einsum('KR,kRxab->kKxab', phase2, #Fourier transform from r2 -> k2
            np.einsum('kr,rRxab->kRxab', phase1.conj(), #Fourier transform from r1 -> k1
            HePhWannier))))) * normFac[...,None,None] #Phonon amplitude factor
        return g, omegaPh, E1, E2, v1, v2
    def returnHwannier():
        return Hwannier
    def returnHeph():
        return HePhWannier
    def returnPwannier():
        return Pwannier
    """
    #Choose a random phonon branch and choose two random electron bands
    n = rand(1:5)
    m = rand(1:5)
    ν = rand(1:3)
    k1 = rand(1, 3)
    k2 = rand(1, 3)
    dir = "../data/electron_phonon_reference/"
    filebase = dir*"wannier"

    export_hwannier(filebase, [12, 12, 12], spin = Val('n'))
    export_heph(filebase, [2, 2, 2], spin = Val('n'))
    export_momentum(filebase, [12, 12, 12], spin = Val('n'))

    Hwannier, cell_map = hwannier(filebase), np.loadtxt(filebase*".map.txt")
    Hephwannier, celleph_map = hephwannier(filebase), np.loadtxt(filebase*".mapeph.txt")

    Pwannier = pwannier(filebase)

    forcematrix, cell_mapph = phonon_force_matrix(joinpath(dir, "totalE"))

    allephs = eph_matrix_elements(Hephwannier, celleph_map, forcematrix, cell_mapph,  Hwannier, cell_map, vec(k1), vec(k2))
    allephspy = py"calcEph"(k1, k2)[1][1, 1, :, :, :]
    @test maximum(abs.(allephs .- 1/eV*allephspy) ./ abs.(allephs) *100) < 5 #Check that maximum difference is less than 5 percent difference. 
    #Check the phonon dispersion. 
    allph = phonon_dispersion(forcematrix, cell_mapph, vec(k1-k2))
    allphpy = py"calcEph"(k1, k2)[2][1, 1, :]*1/eV
    @test maximum((allph .- allphpy)./(allph)*100) < 1 #Check that maximum difference is less than 1 percent
    #Check energies. 
    energypy = vec((py"calcEph"(k1, k2))[3]*1/eV)
    energyjl, _ = wannier_bands(Hwannier, cell_map, vec(k1))
    @test maximum(((energypy .- energyjl)./energyjl)*100) < 1 #Check that energy dispersions are less than 1 percent off

    k1 = rand(10, 3)
    energypy = (py"calcE"(k1))[1]*1/eV
    energyjl, _ = wannier_bands(Hwannier, cell_map, transpose(k1))

    @test energypy ≈ energyjl

    k1 = rand(10, 3)
    k2 = rand(10, 3)

    allephspy = py"calcEph"(k1, k2)[1]*1/eV
    allephs = eph_matrix_elements(Hephwannier, celleph_map, forcematrix, cell_mapph,  Hwannier, cell_map, transpose(k1),transpose(k2))
    @test allephs ≈ allephspy

    @test py"returnHwannier"() ≈ Hwannier
    @test py"returnHeph"() ≈ Hephwannier
    @test py"returnPwannier"() ≈ Pwannier
end