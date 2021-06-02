#=We will include models for the dielectric function of graphene for future reference 
We start first with the density of states of graphene
The functions below are exact expressions for the dynamical polarization of graphene at charge neutrality 
and at finite doping. They are copied verbatim from: B Wunsch et al 2006 New J. Phys. 8 318
=#
function real_neutral(q::Real, w::Real)
    return -1*q^2/(4*(-w^2+36*q^2)^.5)
end

function imag_neutral(q::Real, w::Real)
    return -q^2/(4*(w^2-36*q^2)^.5)
end

function intraband_1a_real(q::Real, w::Real, μ::Real)
    return -2*μ/(36*pi)+1/4*(q^2)/sqrt(abs(36*(q^2)-w^2))
end

function gplus(x::Real)
    @assert x >= 1 "Must supply an argument larger than 1"
    #Note that the log term is a different way of simply writing cosh inverse
    @assert log(x+sqrt(x^2-1)) ≈ acosh(x)
    return x*sqrt(x^2-1)-log(x+sqrt(x^2-1))
end

function gminus(x::Real)
    @assert abs(x) <= 1 "Must give an argument between -1 and 1"
    g = x >= 0 ? x*sqrt(1-x^2)-atan(sqrt(1-x^2)/x) : -pi+x*sqrt(1-x^2)+atan(sqrt(1-x^2)/(-x))
    @assert isapprox(g, x*sqrt(1-x^2) - acos(x), atol=1e-3)
    return g
end

function f(q::Real, w::Real)
    return 1/(4*pi)*q^2/sqrt(abs(w^2-(6*q)^2))
end

function intraband_1b_real(q::Real, w::Real, mu::Real)
    return -2*mu/(6^2*pi)+f(q, w)*(gplus((2*mu+w)/(6*q))-gplus((2*mu-w)/(6*q)))
end

function intraband_1b_imag(q::Real, w::Real, mu::Real)
    return f(q, w)*pi
end

function intraband_1a_imag(q::Real, w::Real, mu::Real)
    return f(q, w)*(gplus((2*mu-w)/(6*q))-gplus((2*mu+w)/(6*q)))
end

function intraband_2a_imag(q::Real, w::Real, mu::Real)
    return -f(q, w)*gplus((2*mu+w)/(6*q))
end

function intraband_2b_imag(q::Real, w::Real, mu::Real)
    return -f(q, w)*gminus((w-2*mu)/(6*q))
end

function intraband_2a_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)-f(q, w)*gminus((w-2*mu)/(6*q))
end

function intraband_2b_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)+f(q, w)*gplus((w+2*mu)/(6*q))
end

function intraband_3a_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)+f(q, w)*(gminus((2*mu+w)/(6*q))-gminus((w-2*mu)/(6*q)))
end

function intraband_3b_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)+f(q, w)*(gplus((2*mu+w)/(6*q))-gplus((w-2*mu)/(6*q)))
end

function intraband_real_total(q::Real, w::Real, mu::Real)
    if w<6*q && w<2*mu-6*q
        return intraband_1a_real(q, w, mu)
    elseif w<-2*mu+6*q
        return intraband_3a_real(q, w, mu)
    elseif w<6*q && w>2*mu-6*q && w>-2*mu+6*q
        return intraband_2a_real(q, w, mu)
    elseif w>6*q && w<2*mu-6*q
        return intraband_1b_real(q, w, mu)
    elseif w>6*q && w>2*mu-6*q && w<2*mu+6*q
        return intraband_2b_real(q, w, mu)
    elseif w>2*mu+6*q
        return intraband_3b_real(q, w, mu)
    else
        return 0
    end
end

function intraband_imag_total(q::Real, w::Real, mu::Real)
    if w<6*q && w<2*mu-6*q
        return intraband_1a_imag(q, w, mu)
    elseif w<-2*mu+6*q #Region 3a
        return 0
    elseif w<6*q && w>2*mu-6*q && w>-2*mu+6*q #Region 2a
        return intraband_2a_imag(q, w, mu)
    elseif w>6*q && w<2*mu-6*q
        return intraband_1b_imag(q, w, mu)
    elseif w>6*q && w>2*mu-6*q && w<2*mu+6*q
        return intraband_2b_imag(q, w, mu)
    elseif w>2*mu+6*q
        return 0
    else
        return 0
    end
end

function graphene_total_polarization(q::Real, w::Real, mu::Real)
    return intraband_real_total(q, w, mu)+ (w<6q ? real_neutral(q, w) : 0 )
end

function graphene_total_impolarization(q::Real, w::Real, mu::Real)
    return intraband_imag_total(q, w, mu) + (w>6q ? imag_neutral(q, w) : 0 )
end

function exact_graphene_epsilon(q::Real, w::Real, mu::Real)
    return 1-e²ϵ/2/q*graphene_total_polarization(q, w, mu) 
end

function exact_graphene_plasmon(q::Real, mu::Real; numevals::Real= 1e6, max_multiple_of_mu::Integer=100, background::Real=1)
    numevals = Int(numevals)
    Epsilons = zeros(numevals)
    for i in 1:numevals
        ω = mu*i/numevals*max_multiple_of_mu
        Epsilons[i] = background-e²ϵ/2/q*graphene_total_polarization(q, ω, mu) 
    end
    return argmin(log.(abs.(Epsilons)))*max_multiple_of_mu/numevals*mu
end

function exact_graphene_plasmonq(ω::Real, mu::Real; numevals::Real=1e6, background::Real=1)
    numevals=Int(numevals)
    logEpsilons=zeros(numevals)
    for i in 1:numevals
        q = mu*i/numevals*200/6
        logEpsilons[i] = log(abs(background-e²ϵ/2/q*(graphene_total_polarization(q, ω, mu))))#+graphene_total_impolarization(q, ω, mu))))
    end
    return argmin(logEpsilons)*200/numevals*mu/6
end

function graphene_plasmon_confinement(λ::Real, μ::Real)
    ω=1.24/λ
    lambdaair=λ*1e-6
    lambdap=2*pi/exact_graphene_plasmonq(ω, μ)*1e-10
    return lambdaair/lambdap
end

"Provides the loss (q2/q1)"
function graphene_plasmon_qloss(λ::Real; μ::Real = 0.135)
    ω = 1.24/λ
    τ = 1.35e-13/ħ ##Tau in units of 1/eV since ω is always given in eV
    q = exact_graphene_plasmonq(ω, μ, background=2.5) #Find the plasmon wavevector
    total_polω =  graphene_total_polarization(q, ω, μ) + 1im*graphene_total_impolarization(q, ω, μ)
    total_pol0 = graphene_total_polarization(q, 0, μ) + 1im*graphene_total_impolarization(q, 0, μ)
    q2_num = graphene_total_impolarization(q, ω, μ) + 1/τ*(graphene_total_polarization(q, ω, μ)-graphene_total_polarization(q, ω-μ/100000, μ))/(μ/100000)+1/(ω*τ)*real(total_polω*(1-total_polω/total_pol0))
    q2_denum = 1/q*graphene_total_polarization(q, ω, μ)-graphene_total_polarization(q, ω, μ)*100000/μ + graphene_total_polarization(q-μ/100000, ω, μ)*100000/μ
    return q/(q2_num/q2_denum)
end

"Provides the group velocity of the graphene plasmon"
function graphene_group_velocity(λ::Real, μ::Real = 0.135)
    ω = 1.24/λ
    q1 = exact_graphene_plasmonq(ω, μ, background=2.5) #Find the plasmon wavevector
    q2 = exact_graphene_plasmonq(ω+μ/30, μ, background=2.5) #Find the plasmon wavevector
    return μ/30/(q2-q1)/(c*ħ)
end

function exact_graphene_landau_damping(q::Real, w::Real, δ::Real, mu::Real)
    RePolω = graphene_total_polarization(q, w, mu) 
    RePolδω = graphene_total_polarization(q, w+δ, mu)  
    ImPol = graphene_total_impolarization(q, w, mu)
    return ImPol*δ/(RePolδω-RePolω)
end

function exact_graphene_landau_damping(q::Real, δ::Real, mu::Real; kwargs...)
    ω = exact_graphene_plasmon(q, mu; kwargs...)
    exact_graphene_landau_damping(q, ω, δ, mu)
end

"""
Provides another method to compute landau damping in graphene, inspired by formalism in the following paper:
Jablan, Marinko, and Darrick E. Chang. "Multiplasmon absorption in graphene." Physical review letters 114.23 (2015): 236801.
"""
function marinko_graphene_landau_damping(q::Real, μ::Real; mesh::Integer = 100, histogram_width::Real=100)
    impolatplas = 0 #Imaginary value of polarization at the plasmon frequency- used to cross check with known values
    loss = 0
    area = 0 
    plasmon = exact_graphene_plasmon(q, μ)
    PlasmonMatrixElement = 4π/137*6.6*3*100*plasmon/(q*4) #4piαhbarc
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        k = i/mesh*μ
        theta = j/mesh*2*π
        kx, ky = k*cos(theta), k*sin(theta)
        kplusq=sqrt((kx+q)^2+ky^2)
        Eupperk, Elowerkplusq = dirac_approximation_upper(k), dirac_approximation_lower(kplusq)
        Eupperkplusq = dirac_approximation_upper(kplusq)
        delta=1
        sameOverlap=1/2*(1+(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5 +delta/100000000))
        mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
        fupperk = heaviside(μ-Eupperk)
        flowerkplusq = heaviside(μ-Elowerkplusq)
        fupperkplusq = heaviside(μ-Eupperkplusq)
        DiffEnergiesUL = Eupperk-Elowerkplusq
        DiffEnergiesUU = Eupperk-Eupperkplusq
        area += (2π/mesh)*k*(μ/mesh)
        if abs(DiffEnergiesUL-plasmon)*histogram_width<0.5 && DiffEnergiesUL>0
            loss += k*(flowerkplusq-fupperk)*mixedOverlap*PlasmonMatrixElement*1/π^2*histogram_width*(μ/mesh)*(2π/mesh)
            impolatplas += -k*(flowerkplusq-fupperk)*mixedOverlap*π/π^2*histogram_width*(μ/mesh)*(2π/mesh)
        end
        if abs(DiffEnergiesUU-plasmon)*histogram_width<0.5 && DiffEnergiesUU>0
            loss += k*(fupperkplusq-fupperk)*sameOverlap*PlasmonMatrixElement*1/π^2*histogram_width*(μ/mesh)*(2π/mesh)
            impolatplas += -k*(fupperkplusq-fupperk)*sameOverlap*π/π^2*histogram_width*(μ/mesh)*(2π/mesh)
        end
    end
    loss *= 2π/ħ
    return loss, impolatplas, graphene_total_impolarization(q, plasmon, μ) 
end

function marinko_graphene_landau_damping_mc(q::Real, μ::Real; mesh::Integer= 100, histogram_width::Integer=100)
    impolatplas = 0 #Imaginary value of polarization at the plasmon frequency- used to cross check with known values
    loss = 0
    plasmon = exact_graphene_plasmon(q, μ, numevals=2000)
    PlasmonMatrixElement = 4π/137*6.6*3*100*plasmon/(q*4) #4piαhbarc
    randcoord = rand(mesh^2, 2)
    for (i, j) in eachrow(randcoord)
        k=i*μ/2
        theta=j*2*π
        kx, ky=k*cos(theta), k*sin(theta)
        kplusq=sqrt((kx+q)^2+ky^2)
        Eupperk, Elowerkplusq = dirac_approximation_upper(k), dirac_approximation_lower(kplusq)
        Eupperkplusq = dirac_approximation_upper(kplusq)
        delta=1
        sameOverlap=1/2*(1+(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5 +delta/100000000))
        mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
        fupperk = heaviside(μ-Eupperk)
        flowerkplusq = heaviside(μ-Elowerkplusq)
        fupperkplusq = heaviside(μ-Eupperkplusq)
        DiffEnergiesUL = Eupperk-Elowerkplusq
        DiffEnergiesUU = Eupperk-Eupperkplusq
        if abs(DiffEnergiesUL-plasmon)*histogram_width<0.5 && DiffEnergiesUL>0
            loss += k*(flowerkplusq)*(1-fupperk)*mixedOverlap*PlasmonMatrixElement*1/π^2*histogram_width*(μ/mesh*0.5)*(2π/mesh)
            impolatplas += -k*(flowerkplusq-fupperk)*mixedOverlap*π/π^2*histogram_width*(μ/(2*mesh))*(2π/mesh)
        end
        if abs(DiffEnergiesUU-plasmon)*histogram_width<0.5 && DiffEnergiesUU>0
            loss += k*(fupperkplusq)*(1-fupperk)*sameOverlap*PlasmonMatrixElement*1/π^2*histogram_width*(μ/mesh*0.5)*(2π/mesh)
            impolatplas += -k*(flowerkplusq-fupperk)*sameOverlap*π/π^2*histogram_width*(μ/(2*mesh))*(2π/mesh)
        end
    end
    loss *= 2π/ħ
    return loss, impolatplas, graphene_total_impolarization(q, plasmon, μ) 
end

function graphene_energy(t::Real, kx::Real, ky::Real)
    t*sqrt(3+2*cos(sqrt(3)*kx*1.42)+4*cos(3/2*ky*1.42)*cos(sqrt(3)*kx/2*1.42))
end

function graphene_energy_normalizedk(t::Real, graphene_lattice::Vector{<:Vector{<:Real}}, k1::Real, k2::Real)
    kx, ky = unnormalize_kvector(graphene_lattice, [k1, k2, 0])
    t*sqrt(3+2*cos(sqrt(3)*kx*1.42)+4*cos(3/2*ky*1.42)*cos(sqrt(3)*kx/2*1.42))
end

function graphene_dos(t::Real, mesh::Real, histogram_width::Real) 
    max_energy=3*abs(t)
    middle_index=round(Int, max_energy*histogram_width)+1
    num_indices=middle_index*2
    GrapheneDOS=zeros(num_indices)
    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]
    K=4*pi/(3*sqrt(3)*1.42);
    N=mesh
    for (i, j) in Tuple.(CartesianIndices(rand(N, N)))
        kxnormal, kynormal=i/N, j/N
        kx, ky = unnormalize_kvector(graphene_lattice, [kxnormal, kynormal, 0])
        Ek=graphene_energy(t, kx, ky) 
        GrapheneDOS[round(Int, histogram_width*Ek)+middle_index]+=(1/N)^2*histogram_width
        GrapheneDOS[-round(Int, histogram_width*Ek)+middle_index]+=(1/N)^2*histogram_width
    end
    return GrapheneDOS
end

function graphene_dos_monte_carlo(t::Real, mesh::Real, histogram_width::Real) 
    max_energy=3*abs(t)
    middle_index=round(Int, max_energy*histogram_width)+1
    num_indices=middle_index*2
    GrapheneDOS=zeros(num_indices)
    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]
    randks = rand(mesh, 2)
    for (rkx, rky) in eachrow(randks)
        kxnormal, kynormal = rkx, rky
        kx, ky = unnormalize_kvector(graphene_lattice, [kxnormal, kynormal, 0])
        Ek=graphene_energy(t, kx, ky)
        GrapheneDOS[round(Int, histogram_width*Ek)+middle_index] += (1/mesh)*histogram_width
        GrapheneDOS[-round(Int, histogram_width*Ek)+middle_index] += (1/mesh)*histogram_width
    end
    return GrapheneDOS
end

function graphene_dos_quad(t::Real, ϵ::Real, δ::Real; kwargs...) 
    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]
    1/π*hcubature(vec->imag(-1/(ϵ-graphene_energy_normalizedk(t, graphene_lattice, vec[1], vec[2])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]
end

function check_graphene_dos_quad(t::Real, δ::Real, npoints::Integer; verbose::Bool=true, kwargs...)     
    quaddos=[]
    for i in 1:npoints
        verbose && println(i)
        ω=i/npoints*abs(t)*3
        push!(quaddos, graphene_dos_quad(t, ω, δ; kwargs...))
    end
    return sum(quaddos*1/npoints*abs(t)*3)
end

"checks that the integrated value of the dos of graphene over all energies gives 2 orbitals per unit cell"
function check_graphene_dos(t::Real, mesh::Real, histogram_width::Real) 
    graphene_dos_array = graphene_dos(t, mesh, histogram_width)
    return sum(graphene_dos_array*1/histogram_width)
end

function dirac_approximation_lower(k::Real)
    return -6*k
end

function dirac_approximation_upper(k::Real)
    return 6*k
end

function lower_band_integrand(k::Real, theta::Real, q::Real, ω::Real, delta::Real)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    return (2*mixedOverlap*(dirac_approximation_lower(k)-dirac_approximation_upper(kplusq))/((dirac_approximation_lower(k)-dirac_approximation_upper(kplusq))^2-(ω+1im*delta)^2))
end

function lower_band_integrand(k::Real, theta::Real, q::Real, ω::Real, μ::Real, delta::Real)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    return (2*mixedOverlap*(dirac_approximation_lowerwself(k, μ)-dirac_approximation_upperwself(kplusq, μ))/((dirac_approximation_lowerwself(k, μ)-dirac_approximation_upperwself(kplusq, μ))^2-(ω+1im*delta)^2))
end

function upper_band_integrand(k::Real, theta::Real, q::Real, ω::Real, delta::Real)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    sameOverlap=1/2*(1+(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5 +delta/100000000))
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    a=2*sameOverlap*(dirac_approximation_upper(k)-dirac_approximation_upper(kplusq))/((dirac_approximation_upper(k)-dirac_approximation_upper(kplusq))^2-(ω+1im*delta)^2)
    b=2*mixedOverlap*(dirac_approximation_upper(k)-dirac_approximation_lower(kplusq))/((dirac_approximation_upper(k)-dirac_approximation_lower(kplusq))^2-(ω+1im*delta)^2)
    return (a+b)
end

function upper_band_integrand(k::Real, theta::Real, q::Real, ω::Real, μ::Real, delta::Real)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    sameOverlap=1/2*(1+(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5 +delta/100000000))
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    eup = dirac_approximation_upperwself(k, μ)
    eupq = dirac_approximation_upperwself(kplusq, μ)
    ednq = dirac_approximation_lowerwself(kplusq, μ)
    a=2*heaviside(μ-real(eup))*sameOverlap*(eup-eupq)/((eup-eupq)^2-(ω+1im*delta)^2)
    b=2*heaviside(μ-real(eup))*mixedOverlap*(eup-ednq)/((eup-ednq)^2-(ω+1im*delta)^2)
    return (a+b)
end

function graphene_conductivity( μ::Real, q::Real, ω::Real; delta::Real=0.01, self::Bool=false, kwargs... )
    if self==false
        A= hcubature(x-> x[1]/(pi^2)*imag(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [μ*3, 2π]; kwargs...)
        B=hcubature(x-> x[1]/(pi^2)*imag(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
        return -4im*ω/q^2*(B[1]+A[1])
    else 
        println(ω)
        A= hcubature( x-> x[1]/(pi^2)*imag(lower_band_integrand(x[1], x[2], q, ω, μ, delta)), [0, 0], [μ*3, 2π]; kwargs...)
        B=hcubature( x-> x[1]/(pi^2)*imag(upper_band_integrand(x[1], x[2], q, ω, μ, delta)), [0, 0], [μ/3, 2π]; kwargs...)
        return -4im*ω/q^2*(B[1]+A[1])
    end
end

function graphene_real_conductivity( μ::Real, q::Real, ω::Real; kwargs... )
    delta = .01
    A= hcubature(x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B = hcubature(x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return 4*ω/q^2*(B[1]+A[1])
end

function graphene_epsilon( μ::Real, q::Real, ω::Real; kwargs... )
    delta=.01
    A = hcubature( x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B = hcubature( x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return 1-e²ϵ/q*(B[1]+A[1])
end

function find_graphene_plasmon(μ::Real, q::Real; nomegas::Integer=3, verbose::Bool=true, kwargs...)
    @info "Numerical calculation of graphene plasmon relation, for exact dispersion use exact_graphene_plasmon"
    epsilon_array=Array{Float64, 1}(undef, nomegas)
    verbose && println("q: ", q)
    for i in 1:nomegas
        ω=i/nomegas*2μ
        verbose && println("ω: ", ω)
        epsilon_array[i]=(log∘abs)(graphene_epsilon( μ, q, ω; kwargs... ))
    end
    return argmin(epsilon_array)/nomegas*2μ
end

function graphene_electron_self_energy(ϵ::Real, μ::Real)
    abs(ϵ-μ)>0.2 ? 0.0183*abs(ϵ-sign(ϵ-μ)*0.2) : 0
end

"""
The matrix elements used in this function are taken from:
Park, Cheol-Hwan, et al. "Velocity renormalization and carrier lifetime in graphene from the electron-phonon interaction." Physical review letters 99.8 (2007): 086804.
"""
function graphene_numerical_self_energy(μ::Real; mesh1::Integer=100, mesh2::Integer=100, histogram_width::Real=100, NQs::Integer=50, 
    verbose::Bool=true)
    g = .035*13.605662285137 # The energy provided in the paper is given in Rydberg
    phononEnergy = 0.2 #Take Optical phonon energy 
    a = 1.42*sqrt(3)
    Abz = brillouin_zone_area([[a, 0, 0], [-a/2, sqrt(3)/2*a, 0], [0, 0, 10]])
    SelfEnergyMat=zeros(NQs)
    for ks in 1:NQs
        verbose && println(ks)
        k = (ks-NQs/2)/NQs*0.8 
        E = k*6
        for (i, j) in Tuple.(CartesianIndices(rand(mesh1, mesh2)))
            q, theta=i/mesh1*1, j/mesh2*(2*π) 
            qx, qy=q*cos(theta), q*sin(theta)
            kplusq=sqrt((k+qx)^2+(qy)^2)
            for band in 1:2
                Energy = (band ==1 ? dirac_approximation_upper(kplusq) : dirac_approximation_lower(kplusq))
                Occupation = heaviside(μ-Energy)
                DiffEnergies1 = Energy+phononEnergy
                DiffEnergies2 = Energy-phononEnergy
                if abs(E-DiffEnergies1)*histogram_width<.5
                    SelfEnergyMat[ks] += (1-Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                end
                if abs(E-DiffEnergies2)*histogram_width<.5
                    SelfEnergyMat[ks] += (Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                end
            end
        end
    end
    return SelfEnergyMat/Abz
end

function graphene_monte_carlo_self_energy(μ::Real; mesh1::Integer=100, mesh2::Integer=100, histogram_width::Real=100, NQs::Integer=50)    
    g = .035*13.605662285137 # The energy provided in the paper is given in Rydberg
    phononEnergy = 0.2 
    SelfEnergyMat=zeros(NQs)
    for ks in 1:NQs
        k=(ks-NQs/2)/NQs*0.4
        E=k*6
        print(ks); flush(stdout)
        random_ks = rand(mesh1)
        random_thetas = rand(mesh2)
        for rks in random_ks
            for thetas in random_thetas
                q, theta=1*rks, 2*π*thetas
                qx, qy=q*cos(theta), q*sin(theta)
                kplusq=sqrt((k+qx)^2+(qy)^2)
                for band in [1, 2]
                    if band==1
                        Energy = dirac_approximation_upper(kplusq)
                    elseif band==2
                        Energy = dirac_approximation_lower(kplusq)
                    end
                    Occupation=heaviside(μ-Energy)
                    DiffEnergies1=Energy+phononEnergy
                    DiffEnergies2=Energy-phononEnergy
                    if  abs(E-DiffEnergies1)*histogram_width<.5
                        SelfEnergyMat[ks]=SelfEnergyMat[ks]+(1-Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                    end
                    if abs(E-DiffEnergies2)*histogram_width<.5
                        SelfEnergyMat[ks]=SelfEnergyMat[ks]+(Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                    end
                end
            end
        end
    end
    return SelfEnergyMat
end

"""
$(TYPEDSIGNATURES)
Returns plasmonic losses in graphene up to first order in the electron-phonon interaction. 

"""
function graphene_second_order_losses(ω::Real; mesh1::Integer=10, mesh2::Integer=200, μ::Real=0.64, δ::Real=0, verbose::Bool=true, wavelength::Bool=true, histogram_width::Real=100)
    wavelength && (ω = 1.24/ω) #Convert from microns to eV if ω is to be interpreted as a wavelength
    verbose && println("Plasmon Frequency is $(ω)")
    #The q below is for Silica on top, vacuum below (2.5 effective dielectric function )
    #q = 2π*2.5*137/(4*pi*ħ*c*μ)*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units, Note that this is in the small wavevector approximation
    #q = 2π*137/(4*pi*ħ*c*μ)*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units, Note that this is in the small wavevector approximation
    q = exact_graphene_plasmonq(ω, μ)
    PlasmonMatrixElement = sqrt(4π/137*6.6*3*100*ω/(q*4)) #4piαhbarc
    #PlasmonMatrixElement /= exact_graphene_epsilon(q, 0, μ) Possible screening to taken into account?
    PhononMatrixElement = 0.4712 # in eV
    Rate = 0 
    A = 5.238760872572826 #Unit Cell Area of Graphene
    for (i, j) in Tuple.(CartesianIndices(rand(mesh1, mesh2)))
        k = i/mesh1*0.5*μ  
        θ = j/mesh2*2π
        ϕi = θ
        ky, kx = k.*sincos(θ)
        for bandi in [-1, 1] #Initial band may be in either lower or upper bands
            ϵi = 6*bandi*k
            fi = heaviside(μ-ϵi)
            isapprox(0, fi) && continue #Only consider filled initial electronic states
            for (l, m) in Tuple.(CartesianIndices(rand(mesh1, mesh2)))
                k2 = l/mesh1*0.5*μ   
                θ2 = m/mesh2*2π     
                ky2, kx2 = k2.*sincos(θ2)
                bandf = 1 #All electronic final states are in upper band
                #Final State Parameters
                ϵf = 6*bandf*sqrt((kx+k2*cos(θ2)+q)^2+(ky+k2*sin(θ2))^2)
                ϕf = atan(real((ky+ky2)/(kx+kx2+q+.000001im)))
                ff = 1-heaviside(μ-ϵf)  
                isapprox(0, ff) && continue #Only consider empty final electronic states

                ϵm = 6*sqrt((kx+kx2)^2+(ky+ky2)^2) #Intermmediate state in upper band, phonon emitted first
                fm = 1-heaviside(μ-ϵm)  
                ϕm = atan(real((ky+ky2)/(kx+kx2+.000001im))) #We take q to lie on x axis WLOG 
                overlapim = 1/2*(1+bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1+bandf*cis(ϕm-ϕf))
                m1 = PlasmonMatrixElement*PhononMatrixElement*(fm)*(overlapim*overlapmf)/(ϵm-ϵi+0.2-1im*δ)

                ϵm = -6*sqrt((kx+kx2)^2+(ky+ky2)^2) #Intermmediate state in lower band, phonon emitted first
                fm = 1-heaviside(μ-ϵm)  
                overlapim = 1/2*(1-bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1-bandf*cis(ϕm-ϕf))
                m2 = PlasmonMatrixElement*PhononMatrixElement*(fm)*(overlapim*overlapmf)/(ϵm-ϵi+0.2-1im*δ)

                ϵm = 6*sqrt((kx+q)^2+ky^2) #Intermmediate state in upper band, plasmon absorbed first
                fm = 1-heaviside(μ-ϵm)  
                ϕm = atan(real((ky)/(kx+q+.000001im))) #We take q to lie on x axis WLOG 
                overlapim = 1/2*(1+bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1+bandf*cis(ϕm-ϕf))
                m3 = PlasmonMatrixElement*PhononMatrixElement*(fm)*(overlapim*overlapmf)/(ϵm-ϵi-ω-1im*δ)

                ϵm = -6*sqrt((kx+q)^2+ky^2) #Intermmediate state in lower band plasmon absorbed first
                fm = 1-heaviside(μ-ϵm)  
                overlapim = 1/2*(1-bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1-bandf*cis(ϕm-ϕf))
                m4 = PlasmonMatrixElement*PhononMatrixElement*(fm)*(overlapim*overlapmf)/(ϵm-ϵi-ω-1im*δ)

                abs(ϵf-ϵi-ω+0.2)*histogram_width<0.5 || continue #Energy conserving delta function
                Rate +=  (4*A/(4*π^2)^2)*k*k2*histogram_width*((2*π*0.5*μ)^2)*1/((mesh1*mesh2)^2)*(abs(m1+m2+m3+m4))^2
            end
        end
    end
    #verbose && println("Landau Damping: ", exact_graphene_landau_damping(q, 0.01, μ)/ħ)
    #verbose && println("Rate: ", Rate)
    return Rate*π/ħ  #+ exact_graphene_landau_damping(q, μ/10, μ, max_multiple_of_mu=100, numevals=1e6)/ħ
end

function graphene_electron_real_self_energy(ϵ::Real, μ::Real, W::Real=8.4)
    pyintegrate.quad(x-> -graphene_electron_self_energy(x, μ)/π, -W, W,  wvar=ϵ, weight="cauchy", limit=10000, epsrel=1e-15, epsabs=1e-15)[1]
end

function dirac_approximation_upperwself(k, μ)
    #6*k+graphene_analytic_real_self_energy(6*k, μ) + 1im*graphene_electron_self_energy(6*k, μ)
    6*k+0.8*ReS(6*k/0.8) + 1im*0.8*ImS(6*k/0.8)
end

function dirac_approximation_lowerwself(k, μ)
   # -6*k+graphene_analytic_real_self_energy(-6*k, μ) + 1im*graphene_electron_self_energy(-6*k, μ)
   -6*k+0.8*ReS(-6*k/0.8) + 1im*0.8*ImS(-6*k/0.8)
end

"""
$(TYPEDSIGNATURES)
Returns the real part of graphene's self energy- corresponding to the band energy shifts due to the electron-phonon interaction at lowest order. 
"""
function graphene_analytic_real_self_energy(ϵ::Real, μ::Real, W::Real=8.4)
    G=0.0183; #Electron-Phonon coupling strength
    w0=0.2; #Phonon frequency
    return G/pi*(w0*log(real(abs((ϵ+.000001im+w0)^2 /(((ϵ+.000001im)-μ)^2-w0^2) ))) - ϵ*log(real(abs(W^2*(ϵ+.000001im-μ+w0)/(((ϵ+.000001im)+w0)^2*(ϵ+.000001im-μ-w0))  )))  );
end

alevitov= 134/sqrt(3); ##Effective nearest neighbor length in angstrom. 134 angstroms is the superlattice size
Klevitov=4*pi/(3*sqrt(3)*alevitov); ##The K point for the superlattice. 
heaviside(x)= x>0 ? 1 : 0

#=
The functions below are to reproduce results from the following paper: 
Intrinsically undamped plasmon modes in narrow electron bands
Cyprian Lewandowski, Leonid Levitov
Proceedings of the National Academy of Sciences Oct 2019, 116 (42) 20869-20874; DOI: 10.1073/pnas.1909069116
=#
function levitov_epsilon(qx::Real, qy::Real, ω::Real; kwargs...)
    q=sqrt(qx^2+qy^2)
    1-e²ϵ*1000/(2*q)*hcubature( x->levitov_integrand(x[1], x[1], qx, qy, ω, .1)*heaviside(limit_up_levitov(x[1])-x[2])*heaviside(-limit_dn_levitov(x[1])+x[2]), [-Klevitov, -Klevitov], [Klevitov, Klevitov]; kwargs...)[1]
end

function limit_dn_levitov(x::Real)
    A1=heaviside(x+Klevitov)*heaviside(-x-Klevitov/2)*(-sqrt(3)*x-4*pi./(3*alevitov));
    A2=heaviside(x+Klevitov/2)*heaviside(-x+Klevitov/2)*(-4*pi/(6*alevitov));
    A3=heaviside(x-Klevitov/2)*heaviside(Klevitov-x)*(sqrt(3)*x-4*pi/(3*alevitov));
    A1+A2+A3;
end

function limit_up_levitov(x::Real)
    B1=heaviside(x+Klevitov)*heaviside(-x-Klevitov/2)*(sqrt(3)*x+4*pi/(3*alevitov));
    B2=heaviside(x+Klevitov/2)*heaviside(-x+Klevitov/2)*(4*pi./(6*alevitov));
    B3=heaviside(x-Klevitov/2)*heaviside(Klevitov-x)*(-sqrt(3)*x+4*pi/(3*alevitov));
    B1+B2+B3;
end

function levitov_energy(kx::Real, ky::Real)
    3.75/3*abs(exp(alevitov*ky*1im)+exp(-(alevitov*kx*sqrt(3)/2+alevitov*ky/2)*1im)+exp((alevitov*kx*sqrt(3)/2-alevitov/2*ky)*1im));
end

function levitov_same_overlap(kx::Real, ky::Real, qx::Real, qy::Real)
    kplusqy=ky+qy;
    kplusqx=kx+qx;
    E1=exp(alevitov*ky*1im)+exp(-(alevitov*kx*sqrt(3)/2+alevitov*ky/2)*1im)+exp((alevitov*kx*sqrt(3)/2-alevitov/2*ky)*1im);
    E2=exp(alevitov*kplusqy*1im)+exp(-(alevitov*kplusqx*sqrt(3)/2+alevitov*kplusqy/2)*1im)+exp((alevitov*kplusqx*sqrt(3)/2-alevitov/2*kplusqy)*1im);
    (1+cos(angle(E1)-angle(E2)))/2;
end

function levitov_mixed_overlap(kx::Real, ky::Real, qx::Real, qy::Real)
    kplusqy=ky+qy;
    kplusqx=kx+qx;
    E1=exp(alevitov*ky*1im)+exp(-(alevitov*kx*sqrt(3)/2+alevitov*ky/2)*1im)+exp((alevitov*kx*sqrt(3)/2-alevitov/2*ky)*1im);
    E2=exp(alevitov*kplusqy*1im)+exp(-(alevitov*kplusqx*sqrt(3)/2+alevitov*kplusqy/2)*1im)+exp((alevitov*kplusqx*sqrt(3)/2-alevitov/2*kplusqy)*1im);
    (1-cos(angle(E1)-angle(E2)))/2
end

function levitov_integrand(kx::Real, ky::Real, qx::Real, qy::Real, w::Real, delta::Real)
    kplusqy=ky+qy;
    kplusqx=kx+qx;
    #=
        The arguments in the exponentials are the "nearest neighbor" distances. 
        The first one is in the y direction with distance alevitov
        The second is at 60 degrees from the negative y axis (sin(30) is 1/2 etc)
        The last one is 30 degrees below the x axis.  
        Note that since the nearest neighbor is in the y direction, the lattice extends in the x direction
        This implies that the reciprocal lattice extends in the y direction. 
    =#
    E1=exp(alevitov*ky*1im)+exp(-(alevitov*kx*sqrt(3)/2+alevitov*ky/2)*1im)+exp((alevitov*kx*sqrt(3)/2-alevitov/2*ky)*1im);
    E2=exp(alevitov*kplusqy*1im)+exp(-(alevitov*kplusqx*sqrt(3)/2+alevitov*kplusqy/2)*1im)+exp((alevitov*kplusqx*sqrt(3)/2-alevitov/2*kplusqy)*1im);
    mixedOverlap=(1-cos(angle(E1)-angle(E2)))/2;
    sameOverlap=(1+cos(angle(E1)-angle(E2)))/2;
    Up1=3.75/3*abs(E1); ## Width of the band is 3.75 meV 
    Up2=3.75/3*abs(E2); 
    a=2*sameOverlap*(Up1-Up2)/((Up1-Up2)^2-(w+1im*delta)^2);
    b=2*mixedOverlap*(Up1+Up2)/((Up1+Up2)^2-(w+1im*delta)^2);
    fullbands1= 1/(12.12*pi^2)*(a+b)*heaviside(1.81-Up1); ##The Fermi energy is at 1.81 eV 
    fullbands2 = 2/(12.12*pi^2)*mixedOverlap*(-Up1-Up2)/((Up1+Up2)^2-(w+1im*delta)^2);

    ###Note that there is an effective background dielectric constant of 12.12 
    ###Note that the prefactors take into account a four fold degeneracy (two layers, two spins)

    fullbands=fullbands1+fullbands2;
    return fullbands
end

function levitov_im_polarization(qx::Real, qy::Real; erange::Real=100, mesh::Integer=100, histogram_width::Integer=100)
    impols = zeros(histogram_width*erange)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(2*mesh+1, 2*mesh+1)))
        x, y = (xmesh-mesh-1)/mesh*Klevitov, (ymesh-mesh-1)/mesh*Klevitov
        (y < limit_up_levitov(x) && y > limit_dn_levitov(x)) || continue
        E1 = levitov_energy(x, y)
        E2 = -levitov_energy(x+qx, y+qy)
        E3 = levitov_energy(x+qx, y+qy)
        #The chemical potential is in the upper band, so we may consider intraband transitions in the upper band
        # and interband transitions between the upper and lower bands
        f1 = heaviside(1.81-E1)
        f2 = heaviside(1.81-E2)
        f3 = heaviside(1.81-E3)
        impols[round(Int, histogram_width*(E1-E2))+1] += (f1-f2)*levitov_mixed_overlap(x, y, qx, qy)*π*(1/π)^2*histogram_width*(Klevitov/mesh)^2
        E1-E3 > 0 || continue
        impols[round(Int, histogram_width*(E1-E3))+1] += (f1-f3)*levitov_same_overlap(x, y, qx, qy)*π*(1/π)^2*histogram_width*(Klevitov/mesh)^2
    end
    return impols
end

function levitov_kramers_kronig_epsilon(qx::Real, qy::Real, ω::Real; kwargs...)
    levitov_impols = levitov_im_polarization(qx, qy; kwargs...)
    histogram_width = 100
    max_energy = 100
    interpolated_ims=interpol.interp1d(0:1/histogram_width:max_energy-1/histogram_width, levitov_impols)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
    q=sqrt(qx^2+qy^2)
    return 12.12-e²ϵ*1000/(2*q)*pyintegrate.quad(cauchy_inner_function, 0, 50, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω)[1]
end

function levitov_kramers_kronig_epsilon(qx::Real, qy::Real, ωs::Vector{<:Real}; kwargs...)
    levitov_impols = levitov_im_polarization(qx, qy; kwargs...)
    histogram_width = 100
    max_energy = 100
    interpolated_ims=interpol.interp1d(0:1/histogram_width:max_energy-1/histogram_width, levitov_impols)
    ErrorAbs=1e-20
    real_epses = Float64[]
    for ω in ωs
        cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
        q=sqrt(qx^2+qy^2)
        push!(real_epses, 12.12-e²ϵ*1000/(2*q)*pyintegrate.quad(cauchy_inner_function, 0, 50, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω)[1])
    end
    return real_epses
end

#=
Next we will examine models for the plasmon modes of bilayer graphene. In this case, the most convenient thing to do 
is to consider the conductivity and solve maxwell's equations. 
=#
"""
$(TYPEDSIGNATURES)
Returns the exact plasmon modes of bilayer graphene (regular not twisted).
Equations used are from:
Gonçalves, Paulo André Dias. Plasmonics and Light–Matter Interactions in Two-Dimensional Materials and in Metal Nanostructures: Classical and Quantum Considerations. Springer Nature, 2020.
"""
function graphene_bilayer_plasmon_modes(q::Real, μ::Real, d::Real; numevals::Real=1e5, max_multiple_of_mu::Integer = 3, background_dielectric::Real=2.5)
    numevals = Int(numevals)
    Diffs=zeros(numevals)
    for i in 1:numevals
        ω = μ*i/numevals*max_multiple_of_mu
        Condoverϵ₀ =  1im*ω/ħ*e²ϵ/q^2*graphene_total_polarization(q, ω, μ)  ##The conductivity divided by ϵ₀
        Diffs[i] = (2*background_dielectric/q+1im*Condoverϵ₀/ω*ħ)^2*exp(2*q*d)-(1im*Condoverϵ₀/ω*ħ)^2
    end
    return log.(abs.(Diffs))
    #return argmin(log.(abs.(Diffs)))*max_multiple_of_mu/num_evals*μ
end

"""
$(TYPEDSIGNATURES)
Numerical calculation of bilayer graphene (non-twisted) plasmonic modes. 

## KEYWORD ARGUMENTS

max_multiple_of_mu : Maximum multiple of the chemical potential to look for plasmon modes

background_dielectric : The constant dielectric function corresponding to the background of the bulk material that surrounds the graphene planes
"""
function find_graphene_bilayer_plasmon_modes(q::Real, μ::Real, d::Real; numevals::Integer = 100, max_multiple_of_mu::Integer = 3, background_dielectric::Real = 2.5, kwargs...)
    delta = 0.01
    Diffs=zeros(numevals)
    for i in 1:numevals
        ω = μ*i/numevals*max_multiple_of_mu
        A = hcubature( x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)[1]
        B = hcubature( x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)[1]
        Condoverϵ₀ =  1im*ω/ħ*e²ϵ/q^2*(A+B)  ##The conductivity divided by ϵ₀
        Diffs[i] = (2*background_dielectric/q+1im*Condoverϵ₀/ω*ħ)^2*exp(2*q*d)-(1im*Condoverϵ₀/ω*ħ)^2
    end
    return log.(abs.(Diffs))
end

function find_graphene_bilayer_plasmon_modes(qs::Vector{<:Real}, μ::Real, d::Real; numevals::Integer = 100, max_multiple_of_mu::Integer = 3, background_dielectric::Real = 2.5, verbose::Bool=true, kwargs...)
    plasmonheatmap = zeros(length(qs), numevals)
    for (idx, q) in enumerate(qs)
        verbose && println(idx)
        modesatq = find_graphene_bilayer_plasmon_modes(q, μ, d; numevals = numevals, max_multiple_of_mu = max_multiple_of_mu, background_dielectric = background_dielectric, kwargs...)
        plasmonheatmap[idx, :] = smooth(modesatq, win_len=20)
    end
    return plasmonheatmap
end