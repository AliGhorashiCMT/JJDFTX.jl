function lattice_line_number(parsed_file::Vector{<:AbstractString}, header::AbstractString)
    linenumber = 0
    for (index, line) in enumerate(parsed_file)
        contains(line, header) || continue
        linenumber = index
        contains(line, header) && break
    end
    return linenumber
end

function loadlattice(outfile::AbstractString)
    parsed_file = readlines(outfile)
    linenumber = lattice_line_number(parsed_file, "R = ")
    row1, row2, row3 = [parse.(Float64, string.(split(parsed_file[linenumber+idx])[2:4])) for idx in 1:3]
    latticearray = hcat(row1, row2, row3)
    return [row for row in eachrow(latticearray)]*bohrtoangstrom
end

"""
$(TYPEDSIGNATURES)
Loads the reciprocal lattice from a JDFTX output file
"""
function loadreciprocallattice(outfile::AbstractString)
    parsed_file = readlines(outfile)
    linenumber = lattice_line_number(parsed_file, "G =")
    row1, row2, row3 = [parse.(Float64, string.(split(parsed_file[linenumber+idx])[2:4])) for idx in 1:3]
    reciprocal_latticearray = hcat(row1, row2, row3)
    return [col for col in eachcol(reciprocal_latticearray)]/bohrtoangstrom
end

function loadcellvolume(outfile::AbstractString)
    Volume = 0 
    for line in readlines(outfile)
        contains(line, "unit cell volume") || continue
        Volume = parse(Float64, string.(split(line))[5])
        contains(line, "unit cell volume") && break
    end
    return Volume*bohrtoangstrom^3
end

function loadcellarea(outfile::AbstractString)
    Volume = loadcellvolume(outfile)
    Zlength = norm(loadlattice(outfile)[3])
    return Volume/Zlength
end

function reciprocal_vectors(lattice_vectors::Vector{<:Vector{<:Real}}) 
    a1, a2, a3 = lattice_vectors
    V = dot(a1, cross(a2, a3))
    return (2π/V)*[cross(a2, a3), cross(a3, a1), cross(a1, a2)]
end

function in_wigner_seitz(lattice_vectors::Vector{<:Vector{<:Real}}, rvec::Vector{<:Real}; n::Integer=2) 
    vec1, vec2, _ = lattice_vectors
    distances_array = Float64[]
    for (i, j) in Tuple.(CartesianIndices(rand(2*n+1, 2*n+1))) ## Iterate from -n to n 
        (i-n-1==0 && j-n-1==0) && continue
        current_vec = vec1*(i-n-1)+vec2*(j-n-1)
        push!(distances_array, euclidean(current_vec, rvec) )
    end
    return (euclidean(rvec, [0, 0, 0]) < minimum(distances_array))
end

function in_brillouin(lattice_vectors::Vector{<:Vector{<:Real}}, kvec::Vector{<:Real}; n::Integer=2) 
    bvectors = reciprocal_vectors(lattice_vectors)
    vec1, vec2, _ = bvectors
    distances_array = Float64[]
    for (i, j) in Tuple.(CartesianIndices(rand(2*n+1, 2*n+1))) ## Iterate from -n to n 
        (i-n-1==0 && j-n-1==0) && continue
        current_vec = vec1*(i-n-1)+vec2*(j-n-1)
        push!(distances_array, euclidean(current_vec, kvec) )
    end
    return (euclidean(kvec, [0, 0, 0]) < minimum(distances_array))
end

function normalize_kvector(lattice_vectors::Vector{<:Vector{<:Real}}, unnormalized_kvector::Vector{<:Real})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    inv(hcat(b1, b2, b3))*unnormalized_kvector
end

function unnormalize_kvector(lattice_vectors::Vector{<:Vector{<:Real}}, normalized_kvector::Vector{<:Real}) 
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    hcat(b1, b2, b3)*normalized_kvector
end

function unit_cell_area(lattice_vectors::Vector{<:Vector{<:Real}}) 
    a1, a2, _ = lattice_vectors
    return norm(cross(a1, a2))
end

function unit_cell_volume(lattice_vectors::Vector{<:Vector{<:Real}}) 
    a1, a2, a3 = lattice_vectors
    return abs(dot(a1, cross(a2, a3)))
end

function brillouin_zone_volume(lattice_vectors::Vector{<:Vector{<:Real}})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    return abs(dot(b1, cross(b2, b3)))
end

function brillouin_zone_area(lattice_vectors::Vector{<:Vector{<:Real}}) 
    b_vectors = reciprocal_vectors(lattice_vectors)
    b_vectors_2d = Vector{Float64}[]
    for b_vector in b_vectors 
        (b_vector[3] ≈ 0) && push!(b_vectors_2d, b_vector)
    end
    b2d_1, b2d_2 = b_vectors_2d 
    return norm(cross(b2d_1, b2d_2))
end


