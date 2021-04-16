"""Returns the eliashberg spectral function. This function is modeled after http://jdftx.org/EphMatrixElements.html
"""
function eliashberg()

#=We have to sum over the brillouin zone twice, and over two electronic band indices and one phonon band index
We have three delta functions. One which enforces the frequency to be equal to the phonon energy, 
=#

#= Units check: 

=#
end

"For use by the Eliashberg spectral function method above"
function dosatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice::Vector{<:Vector{<:Real}}, nbands::Integer, μ::Real; mesh::Integer = 10, histogram_width::Real=3)
    volume = unit_cell_volume(lattice)
    dos = 0 
    for x_mesh in 1:mesh^3
        ϵs = wannier_bands(Hwannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            if abs(μ-ϵ)*histogram_width < 1
                dos = dos + histogram_width*(1/mesh)^3*(1/volume)
            end
        end
    end
    return dos
end

function vFsquaredatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, nbands::Integer, μ::Real; mesh::Integer=10, histogram_width::Real=3)
    vFsquared = 0 
    numintersections = 0 
    for x_mesh in 1:mesh^3
        arandk = rand(3)
        ϵs = wannier_bands(Hwannier, cell_map, arandk, nbands)
        ps = momentum_matrix_elements(Hwannier, cell_map, Pwannier, arandk)
        for (index, ϵ) in enumerate(ϵs)
            if abs(μ-ϵ)*histogram_width < 1
                numintersections +=1
                vFsquared = vFsquared + sum((abs.(ps[:, index, index])).^2)*(bohrtoangstrom/ħ)^2
                ##Note that to stay in keeping with JDFTX conventions, we reconverted to atomic units
            end
        end
    end
    return sqrt(vFsquared/numintersections)
end


function convertdos(dos::Real)
    ##Conventions of this package are that the dos will be in 1/angstrom^3*1/eV units, to convert to jdftx units, 
    ##we must do the following:
    return dos*bohrtoangstrom^3*1/eV
end