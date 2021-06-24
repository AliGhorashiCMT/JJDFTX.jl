@testset "Electron-Phonon Matrix Element Tests" begin
    allattice = [[2, 2, 0], [0, 2, 2], [2, 0, 2]]
    μ = 0.4/eV
    #We test the equivalence of different methods of obtaining the Eliashbert spectral function and the methods therein
    dir = "../data/momentum_matrix_elements"
    wannierfile = joinpath(dir, "Al_wannierbands.txt")
    cellmapfile = joinpath(dir, "Al_cellmap.txt")
    HWannier, cellmap = hwannier(wannierfile, cellmapfile, 5), np.loadtxt(cellmapfile)
    dos = JJDFTX.density_of_states_montecarlo_3d(HWannier, cellmap, 5, mesh=50, energy_range=40, offset=10, histogram_width=5)
    dosmanuallyatmu = dos[round(Int, (μ+10)*5)]/unit_cell_volume(allattice)
    dosatmufunc = JJDFTX.dosatmu(HWannier, cellmap, allattice, 5, μ, mesh=50)
    dosatmulorentzian = JJDFTX.dosatmulorentzian(HWannier, cellmap, allattice, 5, μ, mesh=50, esmearing=0.1)
    dosatmugaussian = JJDFTX.dosatmugaussian(HWannier, cellmap, allattice, 5, μ, mesh=50, esmearing=0.1)
    @test abs(100*(dosmanuallyatmu-dosatmufunc)/dosatmufunc) < 5 #Less than five percent difference
    @test abs(100*(dosatmulorentzian-dosmanuallyatmu)/dosmanuallyatmu) < 5 #Less than five percent 
    @test abs(100*(dosatmugaussian-dosmanuallyatmu)/dosmanuallyatmu) < 5 #Less than five percent 

    PWannier = pwannier(joinpath(dir, "AlP.txt"), cellmapfile)
    @test 0.68 < vFsquaredatmu(HWannier, cellmap, PWannier, 5, μ, mesh=50, histogram_width=10) < 0.73

    #Return fermi k point tests

    validkpoints, subsample = JJDFTX.returnfermikpoint_lorentzian(HWannier, cellmap, 5, μ, esmearing=.01/eV, mesh=3000)
    @test subsample <= 1
    validkpoints, subsample = JJDFTX.returnfermikpoint_gaussian(HWannier, cellmap, 5, μ, esmearing=.01/eV, mesh=3000)
    @test subsample <= 1
    validkpoints, subsample = JJDFTX.returnfermikpoint(HWannier, cellmap, 5, μ, 1,  mesh=3000)
    @test subsample <= 1
    validkpoints2, subsample2 = JJDFTX.returnfermikpoint(HWannier, cellmap, 5, μ, 10,  mesh=3000)
    @test subsample2 <= subsample 
end


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
"""
#Choose a random phonon branch and choose two random electron bands
n = rand(1:5)
m = rand(1:5)
ν = rand(1:3)
#Choose random phonon kvector2 in the Brillouin Zone. 
k1 = rand(1, 3)
k2 = rand(1, 3)
#Now find the electron phonon matrix element
#Compare with our methods
dir = "../data/electron_phonon_reference/"
filebase = dir*"wannier"
HWannier, cellmap = hwannier(filebase*".txt", filebase*".map.txt", 5), np.loadtxt(filebase*".map.txt")
forcematrix, cellmapph = phonon_force_matrix(joinpath(dir, "totalE"))
heph, ceph = write_eph_matrix_elements("../data/electron_phonon_reference/"*"wannier", 3, [2, 2, 2], Val('n'))
allephs = eph_matrix_elements(heph, ceph, forcematrix, cellmapph,  HWannier, cellmap, vec(k1), vec(k2), 5)
allephspy = py"calcEph"(k1, k2)[1][1, 1, :, :, :]
@test maximum(abs.(allephs .- 1/eV*allephspy) ./ abs.(allephs) *100) < 5 #Check that maximum difference is less than 5 percent difference. 

#Check the phonon dispersion. 

allph = phonon_dispersion(forcematrix, cellmapph, vec(k1-k2))
allphpy = py"calcEph"(k1, k2)[2][1, 1, :]*1/eV

@test (allph .- allphpy)./(allph)*100 < 5 #Check that maximum difference is less than 5 percent

#Check energies. 
energypy = vec((py"calcEph"(k1, k2))[3]*1/eV)

energyjl = wannier_bands(HWannier, cellmap, vec(k1), 5)

@test maximum(((energypy .- energyjl)./energyjl)*100) < 5 #Check that energy dispersions are less than 5 percent off

end