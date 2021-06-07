function jdftxcomparewannier(k::Vector{<:Real})
    py"""
    import numpy as np
    from scipy.interpolate import interp1d
    def jdftxcompare(k):
        cellMap = np.loadtxt("../data/graphene_examples/wannier.graphene.in.mlwfCellMap")[:,0:3].astype(int)
        Wwannier = np.fromfile("../data/graphene_examples/wannier.graphene.in.mlwfCellWeights")
        nCells = cellMap.shape[0]
        nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
        Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)
        for line in open('../data/graphene_examples/graphene.out'):
            if line.startswith('kpoint-folding'):
                kfold = np.array([int(tok) for tok in line.split()[1:4]])
        kfoldProd = np.prod(kfold)
        kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
        #--- Read reduced Wannier Hamiltonian and expand it:
        Hreduced = np.fromfile("../data/graphene_examples/wannier.graphene.in.mlwfH").reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)
        iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
        Hwannier = Wwannier * Hreduced[iReduced]
        Hk = np.tensordot(np.exp((2j*np.pi)*np.dot(k, cellMap.T)), Hwannier, axes=1)
        Ek,Vk = np.linalg.eigh(Hk)
        return Ek
    """
    py"jdftxcompare"(k)
end

function jdftxcomparephonon(k::Vector{<:Real})
    py"""
    import numpy as np
    from scipy.interpolate import interp1d
    cellMap = np.loadtxt("../data/graphene_examples/graphene.in.phononCellMap")[:,0:3].astype(int)
    forceMatrix = np.fromfile("../data/graphene_examples/graphene.in.phononOmegaSq", dtype=np.float64)
    nCells = cellMap.shape[0]
    nModes = int(np.sqrt(forceMatrix.shape[0] / nCells))
    forceMatrix = np.reshape(forceMatrix, (nCells,nModes,nModes))
    def jdftxcompareph(k):
        forceMatrixTilde = np.tensordot(np.exp((2j*np.pi)*np.dot(k,cellMap.T)), forceMatrix, axes=1)
        omegaSq, normalModes = np.linalg.eigh(forceMatrixTilde)
        return np.sqrt(omegaSq)*27.2114
    """
    py"jdftxcompareph"(k)
end

@testset "test graphene dispersion" begin
    bands=zeros(8, 30)
    for i in 1:10
        println(i);bands[:, i] = jdftxcomparewannier([0, 0.5*i/10, 0])
    end
    for i in 1:10
        println(i+10);bands[:, i+10] = jdftxcomparewannier([2/3*i/10, 0.5-(0.5+1/3)*i/10, 0])
    end
    for i in 1:10
        println(i+20);bands[:, i+20] = jdftxcomparewannier([2/3-2/3*i/10, -1/3+1/3*i/10, 0])
    end
    @test transpose(bands*1/eV) ≈ JJDFTX.dft_graphene_wannier_dispersion()
end

@testset "test graphene phonon dispersion" begin
    @test dft_graphene_phonon_dispersion([2/3, -1/3, 0]) ≈ jdftxcomparephonon([2/3, -1/3, 0])
    @test dft_graphene_phonon_dispersion([0, -1/2, 0]) ≈ jdftxcomparephonon([0, -1/2, 0])
end