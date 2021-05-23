#using jdftx_to_plot
using Pkg;
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Test, PyCall, JJDFTX

@testset "jdftx_to_plot" begin
    include("wannier_bands_tests.jl")
    include("cellsizes.jl")
    include("dos.jl")
    include("plasmons.jl")
    include("matrix_elements_tests.jl")
    include("analytical_models.jl")
    include("kramers_kronig.jl")
    include("compare_with_jdftx.jl")
end

