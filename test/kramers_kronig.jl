@testset "Simple Interband Model" begin
    ωs = collect(0:1/100:1000)
    real_conductivity = np.heaviside.(ωs.*1/2 .- 1, 0.5)
    ks = [real(kramers_kronig_reverse(ω, real_conductivity, 1000, 0.01)) for ω in 1:1/10:10]
    anals = [-1/2*1/pi*log(((2 .+ ω)./(2 .- ω)).^2) for ω in 1:1/10:10]
    
    for (k, a) in zip(ks, anals)
        (isinf(k) || isinf(a)) && continue
        (isnan(k) || isnan(a)) && continue
        (iszero(k) || iszero(a)) && continue
        @test 100*abs((k-a)/a) < 10 ##Check for less than 10 percent error. 
    end
end