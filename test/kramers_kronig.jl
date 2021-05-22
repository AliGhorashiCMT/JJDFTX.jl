@testset "Simple Interband Model" begin
    ωs = collect(0:1/100:100)
    real_conductivity = np.heaviside.(ωs.*1/2 .- 1, 0.5)
end