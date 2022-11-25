"""
$(TYPEDSIGNATURES)
Loads the lattice from a JDFTX output file

Returns a vector of vectors in angstroms (Note that JDFTx output is in atomic units).

Note that in the case of lattice minimization, the last iteration's lattice sizes are used. 
"""
function loadlattice(outfile::AbstractString)
    linenumber = 0
    parsed_file = readlines(outfile)
    for (index, line) in enumerate(parsed_file)
        contains(line, "R = ") || continue
        linenumber = index
    end
    row1 = parse.(Float64, string.(split(parsed_file[linenumber+1])[2:4]))
    row2 = parse.(Float64, string.(split(parsed_file[linenumber+2])[2:4]))
    row3 = parse.(Float64, string.(split(parsed_file[linenumber+3])[2:4]))
    latticearray = hcat(row1, row2, row3)
    return [row for row in eachrow(latticearray)]*bohrtoangstrom
end

"""
$(TYPEDSIGNATURES)
Loads the reciprocal lattice from a JDFTX output file
"""
function loadreciprocallattice(outfile::AbstractString)
    linenumber = 0
    for (index, line) in enumerate(readlines(outfile))
        line == "G =" || continue
        linenumber = index
    end
    row1 = parse.(Float64, string.(split(readlines(outfile)[linenumber+1])[2:4]))
    row2 = parse.(Float64, string.(split(readlines(outfile)[linenumber+2])[2:4]))
    row3 = parse.(Float64, string.(split(readlines(outfile)[linenumber+3])[2:4]))
    reciprocal_latticearray = hcat(row1, row2, row3)
    return [col for col in eachcol(reciprocal_latticearray)]/bohrtoangstrom
end

"""
$(TYPEDSIGNATURES)
Load the unit cell volume directly from JDFTX output file and convert to angstroms^3
"""
function loadcellvolume(outfile::AbstractString)
    Volume = 0 
    for line in readlines(outfile)
        contains(line, "unit cell volume") || continue
        Volume = parse(Float64, string.(split(line))[5])
        contains(line, "unit cell volume") && break
    end
    return Volume*bohrtoangstrom^3
end

"""
$(TYPEDSIGNATURES)
Load the unit cell area in angstroms^2. Note that the convention is that the third lattice vector is the one in the z direction.

"""
function loadcellarea(outfile::AbstractString)
    Volume = loadcellvolume(outfile)
    Zlength = norm(loadlattice(outfile)[3])
    return Volume/Zlength
end

"""Returns the reciprocal lattice vectors when supplied with three real space vectors

#Examples
```julia-repl
julia> reciprocal_vectors([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
([6.283185307179586, 0.0, 0.0], [0.0, 6.283185307179586, 0.0], [0.0, 0.0, 6.283185307179586])
```
```julia-repl
julia> reciprocal_vectors([[1, 0, 0], [-1/2, √3/2, 0], [0, 0, 1]])
([6.283185307179586, 3.6275987284684357, -0.0], [0.0, 7.255197456936871, 0.0], [0.0, -0.0, 0.6283185307179586])
```
"""
function reciprocal_vectors(lattice_vectors::Vector{<:Vector{<:Real}}) 
    a1, a2, a3 = lattice_vectors
    V = dot(a1, cross(a2, a3))
    return (2π/V)*[cross(a2, a3), cross(a3, a1), cross(a1, a2)]
end

"""
Return whether the vector supplied is within the wigner seitz unit cell.

"""
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


"""Returns the normalized kvector (in the basis of the reciprocal lattice vectors)

```
julia-repl
julia> a = 1.42*√3; K=[4π/(3*a), 0, 0]
julia> graphene_lattice = [[a, 0, 0], [-a/2, a*√3/2, 0], [0, 0, 20]]
julia> normalize_kvector(graphene_lattice, [K, 0, 0])
3-element Array{Float64,1}:
  0.6666666666666665
 -0.33333333333333326
  0.0
```
"""
function normalize_kvector(lattice_vectors::Vector{<:Vector{<:Real}}, unnormalized_kvector::Vector{<:Real})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array = hcat(b1, b2, b3)
    inv(vectors_array)*unnormalized_kvector
end

"""
Returns the wavevector in inverse angstroms when provided the wavevector in the basis of reciprocal lattice vectors


#Examples
```
julia-repl
julia> a = 1.42*√3
julia> graphene_lattice = [[a, 0, 0], [-a/2, a*√3/2, 0], [0, 0, 20]]
julia> unnormalize_kvector(graphene_lattice, [2/3, -1/3, 0])
3-element Array{Float64,1}:
 1.7030979945861202
 0.0
 0.0
 julia> 4π/(3a)
 1.70309799458612
```
"""
function unnormalize_kvector(lattice_vectors::Vector{<:Vector{<:Real}}, normalized_kvector::Vector{<:Real}) 
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array = hcat(b1, b2, b3)
    vectors_array*normalized_kvector
end

"""
Returns the area of the unit cell (required for polarization calculations).

Results are always given in inverse angstroms squared. 
```julia-repl
julia> unit_cell_area([[4, 0, 0], [0, 2, 0], [0, 0, 1]])
8.0
julia> a = 1.42*sqrt(3)
julia> graphene_lattice = lattice([a -a/2 0; 0 a*√3/2 0; 0 0 20])
julia> unit_cell_area(graphene_lattice)
1.467001287260355
```
"""
function unit_cell_area(lattice_vectors::Vector{<:Vector{<:Real}}) 
    a1, a2, _ = lattice_vectors
    return norm(cross(a1, a2))
end


"""
Returns the volume of the unit cell in question. 

#Examples

```julia-repl 
julia> a =1
julia> face_centered_vectors = [[a/2, a/2, 0], [0, a/2, a/2], [a/2, 0, a/2]]
3-element Array{Array{Float64,1},1}:
[0.5, 0.5, 0.0]
[0.0, 0.5, 0.5]
[0.5, 0.0, 0.5]
julia> unit_cell_volume(face_centered_vectors)
0.25
julia> body_centered_vectors = [[a/2, a/2, a/2], [a/2, -a/2, -a/2], [-a/2, -a/2, a/2 ]]
3-element Array{Array{Float64,1},1}:
[0.5, 0.5, 0.5]
[0.5, -0.5, -0.5]
[-0.5, -0.5, 0.5]
julia> unit_cell_volume(body_centered_vectors)
0.5
```

Note that these results comport with our understanding that a body centered cubic has two atoms per conventional unit cell
and a face centered cubic has 4 atoms per conventional unit cell. 

"""
function unit_cell_volume(lattice_vectors::Vector{<:Vector{<:Real}}) 
    a1, a2, a3 = lattice_vectors
    return abs(dot(a1, cross(a2, a3)))
end

function brillouin_zone_volume(lattice_vectors::Vector{<:Vector{<:Real}})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    return abs(dot(b1, cross(b2, b3)))
end


"""Returns the 2d brillouin zone area of the lattice. The assumption is made that the lattice is in the x-y plane
Note that this is equivalent to just dividing 4π^2 by the corresponding unit cell area. 
```julia-repl
julia> a = 1.42*sqrt(3)
julia> graphene_lattice = [[a, 0, 0], [-a/2, a*√3/2, 0], [0, 0, 20]]
julia> brillouin_zone_area(graphene_lattice)
julia> 7.535831194556713
julia> unit_cell_area(graphene_lattice)^-1*(4π^2)
julia> 7.535831194556713
```
"""
function brillouin_zone_area(lattice_vectors::Vector{<:Vector{<:Real}}) 
    b_vectors = reciprocal_vectors(lattice_vectors)
    b_vectors_2d = Vector{Float64}[]
    for b_vector in b_vectors 
        (b_vector[3] ≈ 0) && push!(b_vectors_2d, b_vector)
    end
    b2d_1, b2d_2 = b_vectors_2d 
    return norm(cross(b2d_1, b2d_2))
end


