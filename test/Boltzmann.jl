@testset "Checking Intraband Conductivity of Graphene at Low Doping" begin
   dir = "../data/boltzmann/"
   HWannier, cellmap = hwannier(dir*"wannier.in.txt", dir*"wannier.in.map.txt", 5), np.loadtxt(dir*"wannier.in.map.txt")
   PWannier = pwannier(dir*"p.txt", dir*"wannier.in.map.txt", 5)
   a=1.42
   graphene_lattice=[[a*sqrt(3), 0, 0], [-a*sqrt(3)/2, a*3/2, 0], [0, 0, 10]]
   numerical = drude_conductivity(graphene_lattice, HWannier, cellmap, PWannier, 5,-3.2, mesh=200, histogram_width=3)
   check_against = 4*1/(π) #In terms of the universal conductivity
   @test (100*(numerical-check_against)/numerical) < 10 #Check for less than 10 percent difference 

end

#=
@testset "Checking graphene velocity" begin
   dir = "../data/boltzmann/"
   HWannier, cellmap = hwannier(dir*"wannier-oneband.txt", dir*"wannier-oneband.map.txt", 1), np.loadtxt(dir*"wannier-oneband.map.txt")
   PWannier = pwannier(dir*"p-oneband.txt", dir*"wannier-oneband.map.txt", 1)

   m = 0.5*1e6/(3*1e18)^2
   a=1.42
   graphene_lattice=[[a*sqrt(3), 0, 0], [-a*sqrt(3)/2, a*3/2, 0], [0, 0, 10]]
   #=
   npoints = 200
   velocitiesx = zeros(npoints, npoints)
   velocitiesy = zeros(npoints, npoints)
``
   for (x, y) in Tuple.(CartesianIndices(rand(npoints, npoints)))
      kx, ky, _ = normalize_kvector( graphene_lattice, [(x-npoints/2)/npoints, (y-npoints/2)/npoints, 0]) 
      velocitiesx[x, y] = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, [2/3, -1/3, 0]+[kx, ky, 0]))[1]*ħ/m
      velocitiesy[x, y] = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, [2/3, -1/3, 0]+[kx, ky, 0]))[2]*ħ/m
   end
   =#
   npoints = 1000
   rangetheta = 0:2π/npoints:2π
   velocitiesx = zeros(length(rangetheta))
   velocitiesy = zeros(length(rangetheta))
   k = 1/12
   for (i, θ) in enumerate(rangetheta)
      kx, ky, _ = normalize_kvector( graphene_lattice, [k*cos(θ), k*sin(θ), 0]) 
      velocitiesx[i] = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, [2/3, -1/3, 0]+[kx, ky, 0]))[1]*ħ/m
      velocitiesy[i] = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, [2/3, -1/3, 0]+[kx, ky, 0]))[2]*ħ/m
   end

   (velocitiesy .- sin.(rangetheta))

end

=#