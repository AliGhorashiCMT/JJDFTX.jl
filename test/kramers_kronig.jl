@testset "Simple Interband Model" begin
    ωs = collect(0:1/100:1000)
    real_conductivity = np.heaviside.(ωs.*1/2 .- 1, 0.5)
    ks = [real(kramers_kronig_reverse(ω, real_conductivity, 1000, 0.01)) for ω in 1:1/10:10]
    ksscipy = [real(kramers_kronig_reverse_scipy(ω, real_conductivity, 1000, 0.01, 500)) for ω in 1:1/10:10]
    anals = [-1/2*1/pi*log(((2 .+ ω)./(2 .- ω)).^2) for ω in 1:1/10:10]
    for (k, a) in zip(ks, anals)
        (isinf(k) || isinf(a)) && continue
        (isnan(k) || isnan(a)) && continue
        (iszero(k) || iszero(a)) && continue
        @test 100*abs((k-a)/a) < 10 ##Check for less than 10 percent error. 
    end
    for (k, a) in zip(ksscipy, anals)
        (isinf(k) || isinf(a)) && continue
        (isnan(k) || isnan(a)) && continue
        (iszero(k) || iszero(a)) && continue
        @test 100*abs((k-a)/a) < 10 ##Check for less than 10 percent error. 
    end
    #anals = [-1/2*1/pi*log(((2 .+ ω)./(2 .- ω)).^2) for ω in 0:1/10:10-1/10]
    #replace!(x->isinf(x) ? 0 : x, anals)
    #kramers_kronig(ω, ks, 1000, 1/10) 
end