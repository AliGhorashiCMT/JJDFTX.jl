using Setfield
struct lattice
    rvectors::Array{Float64, 2}
end

function lattice(rvecs::Array{<:Array{<:Real, 1}, 1}) 
    v = Array{Float64, 2}(undef, (3, 3))
    @warn "Note that you created a lattice with dimensions of agnstroms- conversion to Bohr will be automatic"
    for col in 1:3
        v[:, col] = rvecs[col]./bohrtoangstrom
    end
    return lattice(v)
end

struct ionpos
    ionpos::Array{Array{Any, 1}, 1}
end

Base.iterate(ionp::ionpos, state = 1) = state > length(ionp.ionpos) ? nothing : (ionp.ionpos[state], state+1) 
Base.iterate(latt::lattice, state = 1) = state > size(latt.rvectors)[2] ? nothing : (latt.rvectors[:, state], state+1)
Base.length(ionp::ionpos) = length(ionp.ionpos)


function Base.show(io::IO, latt::lattice)
    for (index, l) in enumerate(latt)
        println("Vector $index = $l")
    end
end

function Base.show(io::IO, ions::ionpos)
    for (index, ion) in enumerate(ions)
        println("Ion $index = $ion")
    end
end

struct self_consistent_field
    coulomb_interaction::String
    xc::String
    kpoints::Array{Int, 1}
    lattice::lattice
    ionpos::ionpos
    pseudopotential::String
    charge::T where T<:Real
    spintype::String
    magnetization::R where R<:Real
    smearing::S where S<:Real
    dump::String #Array{String, 1}
end

self_consistent_field(xc, kpoints, lattice, ionpos)=self_consistent_field("Periodic", xc, kpoints, lattice, ionpos,"GBRV/\$ID_pbsesol.uspp", 0, "no-spin", 0, 0.00001, "ElecDensity" )

struct non_self_consistent_field
    scf::self_consistent_field
    kpoints::Array{Int, 1}
end

struct wannier_interpolation
    wannier_centers::Array{Array{Any, 1}, 1}
    saveWFNs::Bool
    phonon::Bool
    phononSupercell::Array{Int64, 1}
    innerWindow::Array{T, 1} where T<:Number
    outerWindow::Array{S, 1} where S<:Number
    wannier_minimize::Int
end

struct phonon
    supercell::Array{Int64, 1}
end

wannier_interpolation(wannier_centers, innerWindow, outerWindow)=wannier_interpolation(wannier_centers, false, false, [0, 0, 0], innerWindow, outerWindow, 10000)

#=
We define functions for multiplying the cells by a given size (for ease of the user)
=#

*(lat::lattice, new_size::Int) = lattice(lat.rvectors*new_size)

function *(scf::self_consistent_field, new_size::Int) 
    new_scf = @set scf.lattice = scf.lattice*new_size
    return new_scf
end

function *(scf::self_consistent_field, cell_mult::Array{Int64, 1})
    new_latt, new_ionpos = make_supercell(scf.lattice, scf.ionpos, cell_mult)
    @set scf.lattice = new_latt
    @set scf.ionpos = new_ionpos
end

function ∘(scf::self_consistent_field, defect_mult::Array{<:Any, 1})
    defect_mult[1]::String
    defect_mult[2]::Array{Int, 1}
    new_latt, new_ionpos = make_defectcell(scf.lattice, scf.ionpos, defect_mult[2], defect_mult[1])

    @set scf.lattice = new_latt
    @set scf.ionpos = new_ionpos

end