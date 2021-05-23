@testset "smooth" begin
    noisysin = sin.(range(0, pi, length=100)) .+ rand(100)./10 .- 1/10
    @test smooth(noisysin, win_len=1) ≈ noisysin
    @test smooth(noisysin, win_len=1000) ≈ noisysin 
    @test smooth(noisysin, win_len=1000, win_method=3) ≈ noisysin 
    @test smooth(noisysin, win_len=1000, win_method=1) ≈ noisysin 
    noisysin = sin.(range(0, pi, length=1000)) .+ rand(1000)./10 .- 1/10
    @test abs((smooth(noisysin, win_len=300) .- sin.(range(0, pi, length=1000)))[500])*100 <= 10 #Less than 10 percent difference from unnoisy data
end