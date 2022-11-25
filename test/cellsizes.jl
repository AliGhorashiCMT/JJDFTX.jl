@testset "cell volumes" begin
    @test unit_cell_volume([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ 1
    @test unit_cell_volume([[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]]) ≈ 1/4
    @test unit_cell_volume([[1/2, 1/2, 1/2], [1/2, -1/2, 1/2], [1/2, 1/2, -1/2]]) ≈ 1/2
end

@testset "cell areas" begin
    @test unit_cell_area([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ 1
    a = 1.42*sqrt(3) 
    @test unit_cell_area([[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]) ≈ a^2*sin(π/3)  
end

@testset "brillouin zone volumes" begin
    @test brillouin_zone_volume([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ (2π)^3
    @test brillouin_zone_volume([[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]]) ≈ 4*(2π)^3
end

@testset "cell brillouin zone areas" begin
    @test brillouin_zone_area([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ 4*π^2
    a = 1.42*sqrt(3) 
    K = 4π/(3a)
    @test brillouin_zone_area([[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]) ≈ 3*sqrt(3)/2*K^2
end

@testset "Normalizing k vectors" begin
    @test normalize_kvector([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [2π, 0, 0]) == [1, 0, 0]
    @test normalize_kvector([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 2π]) == [0, 0, 1]
    @test normalize_kvector([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 2π, 0]) == [0, 1, 0]
    a = 1.42*√3; K = [4π/(3*a), 0, 0];
    graphene_lattice = [[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 20]];
    @test normalize_kvector(graphene_lattice, K) ≈ [2/3, -1/3, 0]
end

@testset "Parsing lattices" begin
    loaded_grlattice = loadlattice((@__DIR__)*"/../data/graphene_examples/graphene.out")
    @test isapprox(loadcellarea((@__DIR__)*"/../data/graphene_examples/graphene.out"), 1.42^2*3*sqrt(3)/2, atol=1e-4)
    @test isapprox(loadcellvolume((@__DIR__)*"/../data/graphene_examples/graphene.out"), 1.42^2*3*sqrt(3)/2*20*.529177, atol=1e-3)
    @test isapprox(sqrt(sum(loadreciprocallattice((@__DIR__)*"/../data/graphene_examples/graphene.out")[2].^2)), 4π/(3*1.42), atol=1e-2)
    @test isapprox(sqrt(sum(loadreciprocallattice((@__DIR__)*"/../data/graphene_examples/graphene.out")[1].^2)), 4π/(3*1.42), atol=1e-2)
end

@testset "Primitive Zones" begin
    @test !in_wigner_seitz([[1, 0, 0], [0, 1, 0], [0, 0, 10]], [0, 0.6, 0])
    @test in_wigner_seitz([[1, 0, 0], [0, 1, 0], [0, 0, 10]], [0, 0, 0])
    @test !in_brillouin([[1, 0, 0], [0, 1, 0], [0, 0, 10]], [2π/2+0.00001, 2π/2-.0001, 0])
    @test in_brillouin([[1, 0, 0], [0, 1, 0], [0, 0, 10]], [2π/2-0.00001, 2π/2-.0001, 0])
end