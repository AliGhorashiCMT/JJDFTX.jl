"Returns the force in eV/angstrom"
function get_force(filename::AbstractString, ion::AbstractString)
    parse.(Float64, split(filter(line -> contains(line, "force $ion "), readlines(filename))[1])[3:5])/eV/bohrtoangstrom
end

"""
Assuming a planar crystal with some out of plane ion, get the distance of the ion from the plane (convention is that 
the plane is at z=0 and that the z direction is perpendicular to the planar crystal)
"""
function get_d(filename::AbstractString, ion::AbstractString)
    lattice_z = loadlattice(filename)[3][3]
    return parse.(Float64, split(filter(line -> contains(line, "ion $ion"), readlines(filename))[1])[5])*lattice_z
 end

 """
 Get the- possibly fixed- initial magnetization.  
 """
 function get_mag(filename::AbstractString)
    magnetization_line = split(first(filter(line->contains(line, "elec-initial-magnetization"), readlines(filename))))
    get_fixed = magnetization_line[3]
    get_fixed == "no" && @warn "Magnetization not fixed- use with care"
    return parse(Float64, magnetization_line[2])
end

"""
Outputs the result of the listEnergy script of jdftx, converted to eV units
"""
function list_energy(filename::AbstractString)
    e = readlines(`listEnergy $filename`)
    e = parse(Float64, first(split(first(e))))
    return e/eV
end

phonon_supercell(filename::AbstractString) = parse.(Int, String.(split(first(filter(line -> contains(line, "phononSupercell"), readlines(filename))))[2:4]))

kpoint_folding(filename::AbstractString) = parse.(Int, String.(split(first(filter(line -> contains(line, "kpoint-folding"), readlines(filename))))[2:4]))

isdone(filename::AbstractString) = contains(last(readlines(filename)), "Done!")

"""
Functionality to check that a particular calculation converged
"""
function is_converged(filename::AbstractString)
    for line in readlines(filename)
        contains(line, "ElecMinimize: Converged") && return true
        contains(line, "SCF: Converged") && return true
    end
    return false
end

"""
$(TYPEDSIGNATURES)

Parse the symmetries outputed in a totalE.sym file. Returns a tuple of rotations and translations (in order).
"""
function parse_symmetries(filebase::AbstractString)
    lines = readlines("$filebase.sym")
    translations = Vector{Vector{Float64}}()
    rotations = Vector{Array{Float64, 2}}()
    rotation = zeros(3, 3)
    translation = zeros(3)
    i = 1
    for line in lines
        parsed_line = parse.(Float64, string.(replace!(split(line),"\"" =>"")))
        if (!isempty(parsed_line) && i < 4)
            rotation[i, :] = parsed_line
            i += 1
        elseif (!isempty(parsed_line) && i == 4)
            translation = round.(parsed_line, digits=3)
        else
            push!(translations, translation)
            push!(rotations, rotation)
            i = 1
            rotation = zeros(3, 3)
            translation = zeros(3)
        end
    end
    return rotations, translations
end


"""
$(TYPEDSIGNATURES)

Find the symmetry opertations that take a kvector to itself modulo a reciprocal lattice vector. 
"""
function find_little_groups(syms::Vector{<:Array{<:Real, 2}}, k::Vector{<:Real}, only_two_d::Bool=true)
    transposed_syms = transpose.(syms)
    only_two_d && (filter!(x->x[3, 3]==1, transposed_syms))
    filter!(x->isapprox(round.(Int, x*k-k ), x*k-k), transposed_syms)
    return transposed_syms
end

"""
$(TYPEDSIGNATURES) 
Returns the symmetry eigenvalue for a particular Kohn Sham wavefunction corresponding to a particular k.

"""
function return_symmetry_eigenvalue(wf::Vector{<:Complex}, ks::Vector{<:Vector{<:Real}}, k::Vector{<:Real},
    iGarr::Vector{<:Vector{<:Vector{<:Real}}}, M::AbstractArray{<:Real, 2}, volume::Real)
    i = findfirst(x->isapprox(x, k, atol=1e-3), ks)
    overlap = 0
    for (idx, kvector) in enumerate(iGarr[i])
        transformed_gk = M*(kvector+k) - k
        transformed_gk = round.(Int, transformed_gk)
        rot_idx = findfirst(x -> isequal(x, transformed_gk), iGarr[i])
        isnothing(rot_idx) && (println("nothing"); continue)
        overlap += conj(wf[idx])*(wf[rot_idx])*volume
    end
    return overlap
end

function return_symmetry_eigenvalue(wf::Vector{<:Complex}, ks::Vector{<:Vector{<:Real}}, k::Vector{<:Real},
    iGarr::Vector{<:Vector{<:Vector{<:Real}}}, M::AbstractArray{<:Real, 2}, v::Vector{<:Real},  volume::Real)
    i = findfirst(x->isapprox(x, k, atol=1e-3), ks)
    overlap = 0
    for (idx, kvector) in enumerate(iGarr[i])
        transformed_gk = M*(kvector+k) 
        rot_idx = findfirst(x -> isequal(x, round.(Int, transformed_gk-k)), iGarr[i])
        isnothing(rot_idx) && (println("nothing"); continue)
        overlap += (wf[idx])*conj(wf[rot_idx])*volume*exp(1im*2*pi*np.dot(v, kvector+k))
    end
    return overlap
end

"""
$(TYPEDSIGNATURES)
"""
function parse_wannier_band_ranges(filebase::AbstractString, spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    filename = 
        if spin == Val('n')
            "$filebase.mlwfBandRanges"
        elseif spin == Val('u')
            "$filebase.mlwfBandRangesUp"
        elseif spin == Val('d')
            "$filebase.mlwfBandRangesDn"
        end
        [parse.(Float64, string.(line)) for line in split.(readlines(filename))]
end

function parse_wannier_centers(filebase::AbstractString, spin::Union{Val{'u'}, Val{'d'}, Val{'n'}}=Val('n'))
    line_nums = Int[]
    wannier_lines = readlines("$filebase.out")
    wannier_centers = Vector{Float64}[]

    for (idx, line) in enumerate(wannier_lines)
        line_num = idx
        contains(line, "Centers in lattice coords:") && push!(line_nums, line_num)
    end

    ((spin isa Val{'u'} || spin isa Val{'d'}) && (length(line_nums) !=2)) && 
    error("Expected output file of spin polarized calculation but found $(length(line_nums)) set of Wannier center(s)")
    
    (spin isa Val{'n'} && (length(line_nums) !=1)) && 
    error("Expected output file of spin polarized calculation but found $(length(line_nums)) set of Wannier center(s)")

    for line_num in line_nums
        for line in wannier_lines[line_num+1:end]
            try
                x, y, z = parse.(Float64, split(line)[2:4])
                push!(wannier_centers, [x, y, z])
            catch
                break
            end
        end
    end
    num_centers = length(wannier_centers)
    wannier_centers = collect(transpose(hcat(wannier_centers...)))
    if spin isa Val{'n'}
        return wannier_centers 
    elseif spin isa Val{'u'} 
        return wannier_centers[1:Int(num_centers/2), :]
    elseif spin isa Val{'d'} 
        return wannier_centers[Int(num_centers/2)+1:num_centers, :]
    end
end