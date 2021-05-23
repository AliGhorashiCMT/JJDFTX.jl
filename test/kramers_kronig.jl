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
    @test isapprox(-1/(2*pi)*log(((2+5)/(2-5))^2), kramers_kronig_reverse(x->heaviside(x/2-1)+heaviside(-x/2-1)-1 , 5, 6), atol=1e-3)
    @test isapprox(kramers_kronig(x->-1/(2*pi)*log(((2+x)/(2-x))^2) , 0.5, 1000), -1, atol=1e-1)
    anals = [-1/2*1/pi*log(((2 .+ ω)./(2 .- ω)).^2) for ω in 0:1/2000:3000]
    replace!(x->isinf(x) ? 0 : x, anals)
    @test isapprox(real(kramers_kronig(0.1, anals, 3000, 2000)), -1, atol=1e-1)
    @test isapprox(kramers_kronig_scipy(0.1, anals, 3000, 2000, 200), 1, atol=1e-1)
    @test isapprox(kramers_kronig_quadgk(0.1, anals, 3000, 2000, 300, δ=0.001), 1, atol=1e-1)
    @test isapprox(kramers_kronig_reverse_quadgk(5, real_conductivity, 1000, 1/100, 10, δ=0.01), -1/(2*pi)*log(((2+5)/(2-5))^2), atol=1e-1)
end