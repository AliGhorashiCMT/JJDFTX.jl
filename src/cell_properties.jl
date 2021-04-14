function cell_vectors(lattice_file::String)
    run(`cat $lattice_file`);
    run(`pwd`)
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
function reciprocal_vectors(lattice_vectors::Array{<:Array{<:Real, 1},1}) 
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    V=dot(a1, cross(a2, a3))
    b1=2π/V*cross(a2, a3)
    b2=2π/V*cross(a3, a1)
    b3=2π/V*cross(a1, a2)
    return b1, b2, b3
end

function reciprocal_vectors(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}}) 
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    V=dot(a1, cross(a2, a3))
    b1=2π/V*cross(a2, a3)
    b2=2π/V*cross(a3, a1)
    b3=2π/V*cross(a1, a2)
    return b1, b2, b3
end

function in_wigner_seitz(lattice_vectors::Array{<:Array{<:Real, 1},1}, rvec::Array{<:Real, 1}) 
    vec1 = lattice_vectors[1]
    vec2 = lattice_vectors[2]
    vec3 = lattice_vectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, rvec) )
            end
        end
    end
    if euclidean(rvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end
end

function in_wigner_seitz(lattice_vectors::lattice, rvec::Array{<:Real, 1}) 
    vec1 = lattice_vectors.rvectors[:, 1]*bohrtoangstrom
    vec2 = lattice_vectors.rvectors[:, 2]*bohrtoangstrom
    vec3 = lattice_vectors.rvectors[:, 3]*bohrtoangstrom
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else 
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, rvec) )
            end
        end
    end
    if euclidean(rvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end
end

function in_brillouin(lattice_vectors::Array{<:Array{<:Real, 1},1}, kvec::Array{<:Real, 1}) 
    bvectors = reciprocal_vectors(lattice_vectors)
    vec1 = bvectors[1]
    vec2 = bvectors[2]
    vec3 = bvectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else 
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, kvec) )
            end
        end
    end
    if euclidean(kvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end
end

function in_brillouin(lattice_vectors::lattice, kvec::Array{<:Real, 1}) 
    lattice_vectors_array = [lattice_vectors.rvectors[:, 1]*bohrtoangstrom,lattice_vectors.rvectors[:, 2]*bohrtoangstrom, lattice_vectors.rvectors[:, 3]*bohrtoangstrom ]
    bvectors = reciprocal_vectors(lattice_vectors_array)
    vec1 = bvectors[1]
    vec2 = bvectors[2]
    vec3 = bvectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else 
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, kvec) )
            end
        end
    end
    if euclidean(kvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end
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
function normalize_kvector(lattice_vectors::Array{<:Array{<:Real, 1},1}, unnormalized_kvector::Array{<:Real, 1})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    inv(vectors_array)*unnormalized_kvector
end

function normalize_kvector(lattice_vectors::Array{<:Array{<:Real, 1},1}, unnormalized_kvector::Tuple{<:Real, <:Real, <:Real})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    inv(vectors_array)*collect(unnormalized_kvector)
end

function normalize_kvector(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}}, unnormalized_kvector::Array{<:Real, 1})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    inv(vectors_array)*unnormalized_kvector
end

function normalize_kvector(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}}, unnormalized_kvector::Tuple{<:Real, <:Real, <:Real})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    inv(vectors_array)*collect(unnormalized_kvector)
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
function unnormalize_kvector(lattice_vectors::Array{<:Array{<:Real, 1},1}, normalized_kvector::Array{<:Real, 1}) 
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    vectors_array*normalized_kvector
end

function unnormalize_kvector(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}}, normalized_kvector::Array{<:Real, 1}) 
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    vectors_array*normalized_kvector
end

function unnormalize_kvector(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}}, normalized_kvector::Tuple{<:Real, <:Real, <:Real}) 
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    vectors_array=Array{Float64,2}(undef, (3, 3))
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3
    vectors_array*collect(normalized_kvector)
end

"""
Returns the area of the unit cell (required for polarization calculations).

Note that if an object of type lattice is passed, the units are assumed to be in Bohr. 

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
function unit_cell_area(lattice_vectors::Array{<:Array{<:Real, 1},1}) 
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    A=sqrt(dot(cross(a1, a2), cross(a1, a2)))
    return A
end

function unit_cell_area(lattice_vectors::lattice)
    a1, a2, a3 = lattice_vectors.rvectors[:, 1]*bohrtoangstrom, lattice_vectors.rvectors[:, 2]*bohrtoangstrom, lattice_vectors.rvectors[:, 3]*bohrtoangstrom
    A=sqrt(dot(cross(a1, a2), cross(a1, a2)))
    return A
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
function unit_cell_volume(lattice_vectors::Array{<:Array{<:Real, 1},1}) 
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    V = abs(dot(a1, cross(a2, a3)))
    return V
end

function unit_cell_volume(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}}) 
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    V = abs(dot(a1, cross(a2, a3)))
    return V
end

"Used as a cross check to make sure the simpler brillouin_zone_volume method is functioning properly"
function brillouin_zone_volume_direct(lattice_vectors::Array{<:Array{<:Real, 1}, 1})
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    VBZ= abs(dot(b1, cross(b2, b3)))
    return VBZ
end

function brillouin_zone_volume_direct(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}} )
    b1, b2, b3 = reciprocal_vectors(lattice_vectors)
    VBZ= abs(dot(b1, cross(b2, b3)))
    return VBZ
end

function brillouin_zone_volume(lattice_vectors::Array{<:Array{<:Real, 1}, 1})
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    V = abs(dot(a1, cross(a2, a3)))
    return (2π)^3/V
end

function brillouin_zone_volume(lattice_vectors::Tuple{Array{<:Real, 1}, Array{<:Real, 1}, Array{<:Real, 1}})
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]
    V = abs(dot(a1, cross(a2, a3)))
    return (2π)^3/V
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
function brillouin_zone_area(lattice_vectors::Array{<:Array{<:Real, 1},1}) 
    b_vectors=reciprocal_vectors(lattice_vectors)
    b_vectors_2d = []
    for b_vector in b_vectors 
        if b_vector[3] ≈ 0
            push!(b_vectors_2d, b_vector)
        end
    end
    b2d_1, b2d_2 = b_vectors_2d 
    return sqrt(sum(cross(b2d_1, b2d_2).^2))
end

function ion_positions(ionpos_file::String)
    run(`cat $ionpos_file`);
    run(`pwd`)
end

function plot_lattice(lattice_file::String)
    ion_position_vectors=String[]
    open(lattice_file, "r") do io
        readline(io)
        ion_position_vectors=readlines(io);
    end
end
