function nonwannier3dimepsilon(filebase::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, numbands::Integer, μ::Real; spin::Integer=1, histogram_width::Real=1)
    momentumfile = "$filebase.momenta"
    eigfile = "$filebase.eigenvals"
    numkpoints = Int(np.shape(np.fromfile(momentumfile, dtype=np.complex))[1]/(3*numbands^2))
    println(numkpoints)
    momenta = np.reshape(np.fromfile(momentumfile, dtype=np.complex), (numkpoints, 3, numbands, numbands))
    momentasquared = sum((abs.(momenta)).^2, dims=2)[:, 1, :, : ] ##Get rid of extra axis now that summation has been performed
    energies = np.reshape(np.fromfile(eigfile), (numkpoints, numbands)) ./ eV
    prefactor = e²ϵ*π*(spin/3)*ħ^2/(mₑ^2) #Factor of three takes into account isotropy. Spin is conventionally taken to be 2
    Epsilons = zeros(100*histogram_width)
    Vol = unit_cell_volume(lattice_vectors)
    for k in 1:numkpoints
        for band1 in 1:numbands
            for band2 in 1:numbands
                energy1, energy2 = energies[k, band1], energies[k, band2]
                ω  = energy2-energy1
                if ω>0
                    pabs = momentasquared[k, band1, band2]*(ħ/bohrtoangstrom)^2
                    f2 = heaviside(μ-energy2)
                    f1 = heaviside(μ-energy1)
                    Epsilons[round(Int, histogram_width*ω)+1] = Epsilons[round(Int, histogram_width*ω)+1] + (f1-f2)*prefactor*(1/Vol)*pabs*histogram_width*(1/ω^2)*1/numkpoints
                end
            end
        end
    end
    return Epsilons
end

#=
The function below will give the imaginary part of the polarization in three dimensions. Note that in order to get epsilon, one need only mutilply by -e^2ϵ/q^2
=#
function nonwannierimpol(filebase::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1},  q::Array{<:Real, 1}, numbands::Integer, μ::Real, ::Val{3}; spin::Integer=1, histogram_width::Real=10)
    momentumfile = "$filebase.momenta"
    eigfile = "$filebase.eigenvals"
    qx, qy, qz = q
    numkpoints = Int(np.shape(np.fromfile(momentumfile, dtype=np.complex))[1]/(3*numbands^2))
    momenta = np.reshape(np.fromfile(momentumfile, dtype=np.complex), (numkpoints, 3, numbands, numbands))
    momentasquared = (abs.(qx*(momenta)[:, 1, :, : ] + qy*(momenta)[:, 2, :, : ] + qz*(momenta)[:, 3, :, : ])).^2##Get rid of extra axis now that summation has been performed
    energies = np.reshape(np.fromfile(eigfile), (numkpoints, numbands)) ./ eV
    Impols = zeros(100*histogram_width)
    V = unit_cell_volume(lattice_vectors)
    for k in 1:numkpoints
        for band1 in 1:numbands
            for band2 in 1:numbands
                band1 == band2 && continue ##Don't consider intraband transitions
                energy1, energy2 = energies[k, band1], energies[k, band2]
                ω  = energy2-energy1
                if ω>0
                    pabs = momentasquared[k, band1, band2]
                    f2 = heaviside(μ-energy2)
                    f1 = heaviside(μ-energy1)
                    overlap = ħ^4/bohrtoangstrom^2*1/(mₑ)^2*1/(ω^2)*momentasquared[k, band1, band2]
                    Impols[round(Int, histogram_width*ω)+1] = Impols[round(Int, histogram_width*ω)+1] + π*(f2-f1)/V*overlap*(1/numkpoints)*histogram_width*spin
                end
            end
        end
    end
    return Impols
end

function nonwannierimpol(filebase::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1},  q::Array{<:Real, 1}, numbands::Integer, μ::Real, ::Val{2}; spin::Integer=1, histogram_width::Real=10)
    momentumfile = "$filebase.momenta"
    eigfile = "$filebase.eigenvals"
    qx, qy, qz = q
    ##Find the number of k points:
    numkpoints = Int(np.shape(np.fromfile(momentumfile, dtype=np.complex))[1]/(3*numbands^2))
    momenta = np.reshape(np.fromfile(momentumfile, dtype=np.complex), (numkpoints, 3, numbands, numbands))
    momentasquared = (abs.(qx*(momenta)[:, 1, :, : ] + qy*(momenta)[:, 2, :, : ] + qz*(momenta)[:, 3, :, : ])).^2##Get rid of extra axis now that summation has been performed
    energies = np.reshape(np.fromfile(eigfile), (numkpoints, numbands)) ./ eV
    Impols = zeros(100*histogram_width)
    V = unit_cell_area(lattice_vectors)
    for k in 1:numkpoints
        for band1 in 1:numbands
            for band2 in 1:numbands
                band1 == band2 && continue ##Don't consider intraband transitions
                energy1, energy2 = energies[k, band1], energies[k, band2]
                ω  = energy2-energy1
                if ω>0
                    pabs = momentasquared[k, band1, band2]
                    f2 = heaviside(μ-energy2)
                    f1 = heaviside(μ-energy1)
                    overlap = ħ^4/bohrtoangstrom^2*1/(mₑ)^2*1/(ω^2)*momentasquared[k, band1, band2]
                    Impols[round(Int, histogram_width*ω)+1] = Impols[round(Int, histogram_width*ω)+1] + π*(f2-f1)/V*overlap*(1/numkpoints)*histogram_width*spin
                end
            end
        end
    end
    return Impols
end