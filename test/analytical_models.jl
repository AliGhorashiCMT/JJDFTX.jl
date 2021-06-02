@testset "graphene conductivity" begin
    @test isapprox(1, imag(JJDFTX.graphene_conductivity(1, 1/50, 3)), atol=1e-1)
    @test isapprox(0, imag(JJDFTX.graphene_conductivity(1, 1/10, 1, atol=1e-10, delta=0.001)), atol=1e-2)
    @test  JJDFTX.exact_graphene_plasmon(1/6, 1) >  JJDFTX.exact_graphene_plasmon(1/12, 1)

    #Below, we will check Kramers Kronig Relations for the imaginary and real polarizations at 0 doping for graphene. 
    q = 1/6
    a = kramers_kronig(ω->JJDFTX.imag_neutral(q, ω), 0.2,  60000, min_energy_integration=1.0000000000001) 
    b = JJDFTX.real_neutral(q, 0.2)
    @test abs((a-b)/a*100) < 1 #Check for less than 1 percent difference

    a = kramers_kronig_reverse(ω->JJDFTX.real_neutral(q, ω), 3, .99999999999)
    b = JJDFTX.imag_neutral(q, 3)
    @test abs((a-b)/a*100)  < 1

    q = 0.05
    a = kramers_kronig(ω->JJDFTX.imag_neutral(q, ω), .06,  6000000, min_energy_integration=.3000001)
    b = JJDFTX.real_neutral(q, 0.06)
    @test abs((a-b)/a*100) < 1 #Check for less than 1 percent difference
    #------------------------------------------------------------------------------------------------------------------------------------#

    #Kramers Kronig Transformations for arbitrary chemical potentials
    for q in 0.1:0.1:1
        a = kramers_kronig(ω->JJDFTX.graphene_total_impolarization(q, ω, 1), 1.3,  60000, min_energy_integration=0)
        b = graphene_total_polarization(q, 1.3, 1)
        @test abs((a-b)/a)*100 < 1 #Less than 1 percent difference
    end

    #Reverse Kramers Kronig Transformations
    for q in 0.1:0.02:3
        a = kramers_kronig_reverse(ω->JJDFTX.graphene_total_polarization(q, ω, 1), 1.3,  60000)
        b = graphene_total_impolarization(q, 1.3, 1)
        (abs(a)>0 && abs(b)>0) && @test (println(q); abs((a-b)/(a))*100 < 1) #Less than 1 percent difference
    end

    @test isapprox(JJDFTX.check_graphene_dos_quad(1, 0.1, 150, maxevals=10000), 1, atol=1e-1)

    a = JJDFTX.exact_graphene_epsilon(1/6, 1.4, 1)
    b = 1-90.5*6*kramers_kronig(ω->JJDFTX.graphene_total_impolarization(1/6, ω, 1), 1.4,  60000, min_energy_integration=0)
    @test abs((a-b)/a*100) < 10 # Less than 10 percent difference


    #Self Energy test 

    selfa = JJDFTX.graphene_electron_self_energy.(range(-2.4, 2.4, length=100), 1)
    self = JJDFTX.graphene_numerical_self_energy(1, NQs=100, mesh1=2000, mesh2=2000,verbose=true)
    selfdif = (self-selfa)
    smalldiff = Integer[]
    for (a, n, d) in zip(selfa, self, selfdif)
        a ≈ 0 && continue
        n ≈ 0 && continue
        d/a*100<5  && push!(smalldiff, 1)#Less than ten percent difference 
    end
    @test sum(smalldiff) > length(self)/2 #Check that for more than half of the cases the difference is small. 
end

@testset "Graphene Self Energy Shift" begin
    numeric = JJDFTX.graphene_electron_real_self_energy.(-3:0.01:3, 1, 5)
    anal = JJDFTX.graphene_analytic_real_self_energy.(-3:0.01:3, 1, 5)
    diffs = (numeric-anal)/numeric
    @test length(findall(x-> abs(x)*100<10, (numeric-anal)./numeric)) > 400
end

@testset "Twisted Bilayer Graphene" begin
    #=ωs = collect(0.5:0.5:20)
    epsilon_levitovs = Array{Float64, 2}(undef, (40, 10))
    Klevitov = JJDFTX.Klevitov
    for (i, k) in enumerate(0:Klevitov/10:Klevitov*9/10)
        println(i); epsilon_levitovs[:, i] = JJDFTX.levitov_kramers_kronig_epsilon(k, 0, ωs )
    end
    =#
    epsilons=zeros(100, 100)
    for (i, k) in enumerate(0:Klevitov/100:Klevitov*99/100)
        println(k)
        for (j, ω) in enumerate(0.2:.2:20)
        epsilons[i, j] = log(abs(JJDFTX.levitov_epsilon(k, 0, ω, maxevals=50)))
        end
    end
    @test 13 <  (0.2:.2:20)[argmin(epsilons[100, :])] < 15 #Check that Plasmon dispersion at K point is between 13 and 15 meV
    @test 13 < (0.2:.2:20)[argmin(log.(abs.(JJDFTX.levitov_kramers_kronig_epsilon(Klevitov, 0, collect(0.2:.2:20)))))] < 15
end

@testset "Graphene Plasmon Scaling" begin
    a = find_graphene_plasmon(1, 1/12, nomegas=20)
    b = find_graphene_plasmon(.1, 1/120, nomegas=20)
    @test isapprox(a, b*10, atol=1e-2)
end

@testset "Monte Carlo DOS sampling" begin
    dosg = JJDFTX.graphene_dos_monte_carlo(1, 450000, 300)
    @test isapprox(sum(dosg)/300, 2, atol=1e-2)
end

@testset "Graphene Landau Damping" begin
    ld = JJDFTX.exact_graphene_landau_damping.(1/120:1/120:1.4/6, 0.001, 1)
    @test isapprox(ld[5], 0, atol=1e-3)
    @test isapprox(ld[10], 0, atol=1e-3)
end

@testset "Graphene Landau Damping With Plasmon Matrix Element" begin
    a, b, c = JJDFTX.marinko_graphene_landau_damping(rand()*2/6, 1, mesh=10000, histogram_width=300)
    if !isapprox(b, 0) && !isapprox(c, 0)
        @test(abs((b-c)/b)*100 < 5 ) #Check for less than five percent error
    end
    a, b, c = JJDFTX.marinko_graphene_landau_damping_mc(rand()*2/6, 1, mesh=5000)
    if !isapprox(b, 0) && !isapprox(c, 0)
        @test(abs((b-c)/b)*100 < 5 ) #Check for less than five percent error
    end
end

@testset "Two Plasmon Absorption in Graphene" begin
    @test isapprox(0, graphenetwoplasmonemission(0.6, 1, mesh=2500, conesize=2, δ=0, histogram_width=10), atol=1e-5) 
end

@testset "Finding Plasmon Wavevector in Graphene For Given Excitation Frequency" begin
    μ = rand(0.1:0.1:1)
    ω = rand((0.1*μ):(0.1*μ):(2*μ))
    q = JJDFTX.exact_graphene_plasmonq( ω, μ, numevals=1e7)
    epsilon = JJDFTX.exact_graphene_epsilon(q, ω, μ)
    @test isapprox(epsilon, 0, atol=1e-3 ) #Check that q corresponds to the zero of the logitudinal dielectric function 
end

@testset "Phonon assisted plasmon losses" begin
    #Check that plasmon frequencies smaller than the optical phonon frequency do not have first order (one phonon 1 electron-hole pair) losses
    #Note that this must be true due to energy conservation. The electron hole pair has positive energy, the phonon has positive energy (being emitted in final state), so the plasmon must have 
    #energy more than 0.2 eV 
    μ = rand(0.1:0.1:1)
    @test isapprox(0, JJDFTX.graphene_second_order_losses.(1.24/0.18, mesh1=30, mesh2=30, μ=μ,  δ=0.01, histogram_width=20)*ħ, atol=1e-4)
end

@testset "Bilayer Graphene Modes" begin
    bilayermodes = find_graphene_bilayer_plasmon_modes(0.4/6*1/2, 0.4, 50, background_dielectric=1, numevals=30, maxevals=100000)
    @test 1 < argmin(bilayermodes)*3/30 < 1.5 #Compare to known values for plasmon dispersion 
    analyticbilayermodes = graphene_bilayer_plasmon_modes(0.4*1/6*1/2, 0.4, 50, background_dielectric=1, numevals=30)
    @test 1 < argmin(analyticbilayermodes)*3/30 < 1.5 
end