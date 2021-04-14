
"""
Returns a supercell with the designated multiplicity from a smaller unit cell 

```
julia-repl
julia> simple_cubic = lattice([1 0 0; 0 1 0; 0 0 1])
julia> small_ionpos = ionpos([["ion", "C", 0, 0, 0, 1], ["ion", "C", 0.5, 0.5, 0.5, 1]])
julia> supercell_ionpos = make_supercell(simple_cubic, small_ionpos, [2, 2, 2])[2]
julia>  for ion in supercell_ionpos.ionpos 
            println(ion)
        end
Any["ion", "C", 0.0, 0.0, 0.0, 1.0]
Any["ion", "C", 0.0, 0.0, 0.5, 1.0]
Any["ion", "C", 0.0, 0.5, 0.0, 1.0]
Any["ion", "C", 0.0, 0.5, 0.5, 1.0]
Any["ion", "C", 0.5, 0.0, 0.0, 1.0]
Any["ion", "C", 0.5, 0.0, 0.5, 1.0]
Any["ion", "C", 0.5, 0.5, 0.0, 1.0]
Any["ion", "C", 0.5, 0.5, 0.5, 1.0]
Any["ion", "C", 0.25, 0.25, 0.25, 1.0]
Any["ion", "C", 0.25, 0.25, 0.75, 1.0]
Any["ion", "C", 0.25, 0.75, 0.25, 1.0]
Any["ion", "C", 0.25, 0.75, 0.75, 1.0]
Any["ion", "C", 0.75, 0.25, 0.25, 1.0]
Any["ion", "C", 0.75, 0.25, 0.75, 1.0]
Any["ion", "C", 0.75, 0.75, 0.25, 1.0]
Any["ion", "C", 0.75, 0.75, 0.75, 1.0]
```
"""
function make_supercell(small_lattice::lattice, small_ionpos::ionpos, cell_mult::Array{Int, 1})
    mult1, mult2, mult3 = cell_mult
    supercell_ionpos = []
    for ion_positions in small_ionpos.ionpos
        new_base_position = ion_positions[3:end].*(1/mult1, 1/mult2, 1/mult3, 1)
        label_1 = ion_positions[1]
        label_2 = ion_positions[2]
        for i in 0:mult1-1
            for j in 0:mult2-1
                for k in 0:mult3-1
                    new_pos=[] 
                    push!(new_pos, label_1, label_2)
                    append!(new_pos, new_base_position +[i/mult1, j/mult2, k/mult3, 0] )
                    push!(supercell_ionpos, new_pos)
                end
            end
        end
    end
    return lattice([small_lattice.rvectors[:, 1]*cell_mult[1] small_lattice.rvectors[:, 2]*cell_mult[2] small_lattice.rvectors[:, 3]*cell_mult[3]]), ionpos(supercell_ionpos) 
end

"""
Returns a defect cell with the designated multiplicity from a smaller unit cell and a designated defect atom

```
julia-repl
julia> simple_cubic = lattice([1 0 0; 0 1 0; 0 0 1])
julia> small_ionpos = ionpos([["ion", "C", 0, 0, 0, 1], ["ion", "C", 0.5, 0.5, 0.5, 1]], "B")
julia> defectcell_ionpos = make_supercell(simple_cubic, small_ionpos, [2, 2, 2])[2]
julia>  for ion in defectcell_ionpos.ionpos
            println(ion)
        end
Any["ion", "B", 0.0, 0.0, 0.0, 1.0]
Any["ion", "C", 0.0, 0.0, 0.5, 1.0]
Any["ion", "C", 0.0, 0.5, 0.0, 1.0]
Any["ion", "C", 0.0, 0.5, 0.5, 1.0]
Any["ion", "C", 0.5, 0.0, 0.0, 1.0]
Any["ion", "C", 0.5, 0.0, 0.5, 1.0]
Any["ion", "C", 0.5, 0.5, 0.0, 1.0]
Any["ion", "C", 0.5, 0.5, 0.5, 1.0]
Any["ion", "C", 0.25, 0.25, 0.25, 1.0]
Any["ion", "C", 0.25, 0.25, 0.75, 1.0]
Any["ion", "C", 0.25, 0.75, 0.25, 1.0]
Any["ion", "C", 0.25, 0.75, 0.75, 1.0]
Any["ion", "C", 0.75, 0.25, 0.25, 1.0]
Any["ion", "C", 0.75, 0.25, 0.75, 1.0]
Any["ion", "C", 0.75, 0.75, 0.25, 1.0]
Any["ion", "C", 0.75, 0.75, 0.75, 1.0]
```
"""
function make_defectcell(small_lattice::lattice, small_ionpos::ionpos, cell_mult::Array{Int, 1}, defect_atom::String)
    defect_lattice, supercell_ionpos =  make_supercell(small_lattice, small_ionpos, cell_mult)
    supercell_ionpos.ionpos[1][2] = defect_atom
    return defect_lattice, supercell_ionpos
end

function make_defectcells(small_lattice::lattice, small_ionpos::ionpos, cell_mults::Array{Array{<:Integer, 1}}, defect_atom::String)
    defect_lattices = Vector{lattice}
    defect_ionposes = Vector{ionpos}
    for cell_mult in cell_mults
        defect_lattice, supercell_ionpos = make_supercell(small_lattice, small_ionpos, cell_mult)
        supercell_ionpos.ionpos[1][2] = defect_atom
        push!(defect_lattices, defect_lattice)
        push!(defect_ionposes, supercell_ionpos)
    end
    return defect_lattices, defect_ionposes
end

function make_bilayer(monolayer_lattice::lattice, monolayer_ionpos::ionpos, distance::Real, translation::Array{<:Real, 1})
    new_ionpos = []
    for (index, ion) in enumerate(monolayer_ionpos.ionpos)
        push!(new_ionpos, ion)
        second_ion = [ion[1], ion[2], ion[3]+translation[1], ion[4]+translation[2], ion[5]+distance, ion[6]]
        push!(new_ionpos, second_ion)
    end
    return monolayer_lattice, ionpos(new_ionpos)
end