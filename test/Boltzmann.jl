@testset "Checking Boltzmann Transport Time of Bulk Silver" begin
   dir = "../data/Silver/"
   Hwannier, cell_map = hwannier("$(dir)wannier"), np.loadtxt("$(dir)wannier.map.txt");
   force_matrix, cellph_map = phonon_force_matrix("$(dir)totalE");
   Pwannier = pwannier("$(dir)wannier");
   Heph, celleph_map = hephwannier("$(dir)wannier"), np.loadtxt("$(dir)wannier.mapeph.txt");
   lattice_vectors = loadlattice("$(dir)wannier.out");

   subsampling = returnfermikpoint(Hwannier, cell_map, 13, Val(3), num_blocks=100, mesh=10, histogram_width=1)
   
   dosmu = JJDFTX.dosatmu(Hwannier, cell_map, 13, Val(3), mesh=10, histogram_width=1, num_blocks=50)
   
   tau =  τ(Hwannier, cell_map, Pwannier, force_matrix, cellph_map,
      Heph, celleph_map, collect(0.01:0.01:1), 13, Val(:histogram), Val(3); 
      histogram_width=1, supplysampling=subsampling, num_blocks=5, supplydos=dosmu, mesh=64, fracroom=1)

   @test 20 < tau[1] < 40
   
end


@testset "Checking Intraband Conductivity of Graphene at Low Doping" begin
   dir = "../data/boltzmann/"

   Hwannier, cell_map = hwannier("$dir/wannier"), np.loadtxt("$dir/wannier.map.txt");
   Pwannier = pwannier("$dir/wannier");   
   lattice_vectors = loadlattice("$dir/graphene.out");
   x, y = drude_conductivity(lattice_vectors, Hwannier, cell_map, Pwannier, mesh=32, num_blocks=100,
   histogram_width=20, degeneracy=2)   
   
   numerical = y[argmin(abs.(x .+ 4.225963500 .- 0.5)), 1, 1] 

   check_against = 2*1/(π) #In terms of the universal conductivity

   @test (100*(numerical-check_against)/numerical) < 10 #Check for less than 10 percent difference 

end

@testset "Checking Interband Conductivity of Graphene at Low Doping" begin
   dir = "../data/boltzmann/"
   dir = "../data/boltzmann/"

   Hwannier, cell_map = hwannier("$dir/wannier"), np.loadtxt("$dir/wannier.map.txt");
   Pwannier = pwannier("$dir/wannier");   
   lattice_vectors = loadlattice("$dir/graphene.out");

   sigma1 = interbandsigma(lattice_vectors, Hwannier, cell_map, Pwannier, -4.22596350037659+0.5, mesh=32,
      num_blocks=500, degeneracy=2, histogram_width=100, energy_range=3)

   @test 0.3 < real((sigma1[1, 1, :]+sigma1[2, 2, :])/8)[end] < 0.5

end

