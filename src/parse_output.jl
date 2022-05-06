# All utilities for parsing jdftx output

"""
Returns the force in eV/angstrom
"""
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
Outputs the result of the listEnergy script of jdftx
"""
function list_energy(filename::AbstractString)
    e = readlines(`listEnergy $filename`)
    e = parse(Float64, first(split(first(e))))
    return e/eV
end

"""
Functionality to check that a particular jdftx calculation finished successfully
"""
function isdone(filename::AbstractString)
    return contains(last(readlines(filename)), "Done!")
end

"""
Functionality to check that a particular calculation converged
"""
function is_converged(filename::AbstractString)
    for line in readlines(filename)
        contains(line, "ElecMinimize: Converged") && return true
        contains(line, "SCF: Converged") && return true
    end
end

function parse_symmetries(filebase::AbstractString)
    lines = readlines("$filebase.sym")
    syms = Vector{Array{Float64, 2}}()
    sym = zeros(3, 3)
    i = 1
    for line in lines
        parsed_line = parse.(Float64, string.(replace!(split(line),"\"" =>"")))
        if !isempty(parsed_line)
            i == 4 && continue
            sym[i, :] = parsed_line
            i += 1
        else
            push!(syms, sym)
            i = 1
            sym = zeros(3, 3)
        end
    end
    return syms
end

function find_little_groups(syms::Vector{<:Array{<:Real, 2}}, k::Vector{<:Real}, only_two_d::Bool=true)
    transposed_syms = transpose.(syms)
    only_two_d && (filter!(x->x[3, 3]==1, transposed_syms))
    filter!(x->isapprox(round.(Int, x*k-k ), x*k-k), transposed_syms)
    return transposed_syms
end

function return_symmetry_eigenvalue(wf::Vector{<:Complex}, ks::Vector{<:Vector{<:Real}}, k::Vector{<:Real},
    iGarr::Vector{<:Vector{<:Vector{<:Real}}}, M::AbstractArray{<:Real, 2}, volume::Real)

    i = findfirst(x->isapprox(x, k, atol=1e-3), ks)
    println(i)
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