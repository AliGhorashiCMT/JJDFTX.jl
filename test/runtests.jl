using Pkg;
using Test, PyCall, JJDFTX

@testset "jdftx_to_plot" begin
    include("matrix_elements_tests.jl")
    include("cellsizes.jl")
    include("smooth.jl")
    include("kramers_kronig.jl")
    include("dos.jl")
    
    #include("plasmons.jl")
    
    #include("analytical_models.jl")
    #include("compare_with_jdftx.jl")
    
    
    #include("heat_capacity.jl")
    #include("density_plots.jl")
    
    #include("Boltzmann.jl")
    
    #include("exporting_hamiltonians.jl")
end


