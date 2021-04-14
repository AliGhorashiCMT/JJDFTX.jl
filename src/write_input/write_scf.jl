function write_lattice(lattice_vectors::lattice, filename::String)
    open(filename, create=true, write=true) do io
        write(io, "lattice \\ \n")
        for lattice_row in eachrow(lattice_vectors.rvectors)
            for coord in lattice_row
                write(io, string(coord), " ")
            end
            write(io, "\\ \n" )
        end
    end
end

function write_lattices(lattice_vector_set::Vector{lattice}, filenames::Vector{<:String})
    #for (index, lattice_vector) in enumerate(lattice_vector_set)
    for (lattice_vector, filename) in zip(lattice_vector_set, filenames)
        open(filename, create=true, write=true) do io
            write(io, "lattice \\ \n")
            for lattice_row in eachrow(lattice_vectors.rvectors)
                for coord in lattice_row
                    write(io, string(coord), " ")
                end
                write(io, "\\ \n" )
            end
        end
    end
end

function write_ionpos(ions::ionpos, filename::String)
    open(filename, create=true, write=true) do io
        for ion in ions.ionpos
            for coord in ion
                write(io, string(coord))
                write(io, "  ")
            end
            write(io, "\n")
        end
    end
end

function write_ionposes(ions_set::Array{ionpos}, filenames::Vector{<:String})
    #for (index, ions) in enumerate(ions)
    for (ions, filename) in zip(ions_set, filenames)
        open(filename, create=true, write=true) do io
            for ion in ions.ionpos
                for coord in ion
                    write(io, string(coord))
                    write(io, "  ")
                end
                write(io, "\n")
            end
        end
    end
end

function write_scf(scf::self_consistent_field, filename::String, ionpos_filename::String, lattice_filename::String)
    open(filename, create=true, write=true, append=false) do io
        write(io, "coulomb-interaction $(scf.coulomb_interaction) \n" )
        write(io, "include  $(ionpos_filename) \n")
        write(io, "include  $(lattice_filename) \n")
        write(io, "ion-species $(scf.pseudopotential)\n")
        write(io, "elec-cutoff 20 100\n")
        write(io, "elec-initial-charge $(scf.charge)\n")
        if scf.spintype != "no-spin"
            write(io, "elec-initial-magnetization $(scf.magnetization) no \n")
        end
        write(io, "spintype $(scf.spintype)\n")
        write(io, "electronic-SCF\n")
        write(io, "dump-name $(string(filename[1:end-3], ".", "\$", "VAR"))\n")
        write(io, "dump End $(scf.dump)\n")
        write(io, "kpoint-folding $(scf.kpoints[1])  $(scf.kpoints[2])  $(scf.kpoints[3])\n")
        write(io, "elec-smearing Fermi $(scf.smearing)\n")
        write(io, "elec-ex-corr $(scf.xc)", "\n")
    end
end

function write_scf(scf::self_consistent_field, filebase::String)
    open("$(filebase).in", create=true, write=true, append=false) do io
        write(io, "coulomb-interaction $(scf.coulomb_interaction) \n" )
        write(io, "include  $(filebase).ionpos \n")
        write(io, "include  $(filebase).lattice \n")
        write(io, "ion-species $(scf.pseudopotential)\n")
        write(io, "elec-cutoff 20 100\n")
        write(io, "elec-initial-charge $(scf.charge)\n")
        if scf.spintype != "no-spin"
            write(io, "elec-initial-magnetization $(scf.magnetization) no \n")
        end
        write(io, "spintype $(scf.spintype)\n")
        write(io, "electronic-SCF\n")
        write(io, "dump-name $(string(filebase, ".", "\$", "VAR"))\n")
        write(io, "dump End $(scf.dump)\n")
        write(io, "kpoint-folding $(scf.kpoints[1])  $(scf.kpoints[2])  $(scf.kpoints[3])\n")
        write(io, "elec-smearing Fermi $(scf.smearing)\n")
        write(io, "elec-ex-corr $(scf.xc)", "\n")
    end
end

function write_scfs(scfs::Vector{self_consistent_field}, filebases::Vector{<:String})
    #for (index, filebase) in enumerate(fielbases)
    for (filebase, scf) in zip(filebases, scfs)
        open("$(filebase).in", create=true, write=true, append=false) do io
            write(io, "coulomb-interaction $(scf.coulomb_interaction) \n" )
            write(io, "include  $(filebase).ionpos \n")
            write(io, "include  $(filebase).lattice \n")
            write(io, "ion-species $(scf.pseudopotential)\n")
            write(io, "elec-cutoff 20 100\n")
            write(io, "elec-initial-charge $(scf.charge)\n")
            if scf.spintype != "no-spin"
                write(io, "elec-initial-magnetization $(scf.magnetization) no \n")
            end
            write(io, "spintype $(scf.spintype)\n")
            write(io, "electronic-SCF\n")
            write(io, "dump-name $(string(filebase, ".", "\$", "VAR"))\n")
            write(io, "dump End $(scf.dump)\n")
            write(io, "kpoint-folding $(scf.kpoints[1])  $(scf.kpoints[2])  $(scf.kpoints[3])\n")
            write(io, "elec-smearing Fermi $(scf.smearing)\n")
            write(io, "elec-ex-corr $(scf.xc)", "\n")
        end
    end
end

function write_kpoints(kvec_coords::Vector{<:Vector{<:Real}}, kvec_labels::Vector{<:AbstractString}, spacing::Real)
    total_kvecs = Vector{Vector{Any}}()
    for (index, coord) in enumerate(kvec_coords)
        push!(total_kvecs, ["kpoint", coord..., kvec_labels[index]])
    end
    open("bandstruct.kpoints.in", "w") do io
            writedlm(io, total_kvecs); write(io, " \n ")
    end;
    run(`bandstructKpoints bandstruct.kpoints.in $(spacing) bandstruct`) 
    rm("bandstruct.kpoints.in")
    rm("bandstruct.plot")
end

function write_randkpoints(numkpoints::Integer, ::Val{3})
    total_kvecs = Vector{Vector{Any}}()
    kpoints = []
    open("bandstruct.randkpoints", "w") do io
        write(io, "kpoint-folding 1 1 1 \n\n\n")
        for i in 1:numkpoints
            kpoint = ["kpoint", rand(3)..., 1/numkpoints]
            push!(kpoints, kpoint)
        end
        writedlm(io, kpoints); write(io, " \n ")
    end;
end

function write_randkpoints(numkpoints::Integer, ::Val{2})
    total_kvecs = Vector{Vector{Any}}()
    kpoints = []
    open("bandstruct.randkpoints", "w") do io
        write(io, "kpoint-folding 1 1 1 \n\n\n")
        for i in 1:numkpoints
            kpoint = ["kpoint", rand(2)..., 0, 1/numkpoints]
            push!(kpoints, kpoint)
        end
        writedlm(io, kpoints); write(io, " \n ")
    end;
end

function write_nscf(scf::self_consistent_field, filename::String, scf_filename::String, ionpos_filename::String, lattice_filename::String, kpoints::String)
    open(filename, create=true, write=true, append=false) do io
        write(io, "coulomb-interaction $(scf.coulomb_interaction) \n" )
        write(io, "include  $(ionpos_filename) \n")
        write(io, "include  $(lattice_filename) \n")
        write(io, "ion-species $(scf.pseudopotential)\n")
        write(io, "elec-cutoff 20 100\n")
        write(io, "elec-initial-charge $(scf.charge)\n")
        if scf.spintype != "no-spin"
            write(io, "elec-initial-magnetization $(scf.magnetization) no \n")
        end
        write(io, "spintype $(scf.spintype)\n")
        write(io, "dump-name $(string(filename[1:end-3], ".", "\$", "VAR"))\n")
        write(io, "dump End $(scf.dump)\n")
        write(io, "include $(kpoints) \n")
        write(io, "fix-electron-density $(scf_filename).\$VAR \n")
        write(io, "elec-smearing Fermi $(scf.smearing)\n")
        write(io, "elec-ex-corr $(scf.xc)", "\n")
    end
end

function write_nscf(scf::self_consistent_field, filebase::String, kpoints::String)
    open("$(filebase).nscf.in", create=true, write=true, append=false) do io
        write(io, "coulomb-interaction $(scf.coulomb_interaction) \n" )
        write(io, "include  $(filebase).ionpos \n")
        write(io, "include  $(filebase).lattice \n")
        write(io, "ion-species $(scf.pseudopotential)\n")
        write(io, "elec-cutoff 20 100\n")
        write(io, "elec-initial-charge $(scf.charge)\n")
        if scf.spintype != "no-spin"
            write(io, "elec-initial-magnetization $(scf.magnetization) no \n")
        end
        write(io, "spintype $(scf.spintype)\n")
        write(io, "dump-name $(string(filebase, ".", "\$", "VAR"))\n")
        write(io, "dump End $(scf.dump)\n")
        write(io, "include $(kpoints) \n")
        write(io, "fix-electron-density $(filebase).\$VAR \n")
        write(io, "elec-smearing Fermi $(scf.smearing)\n")
        write(io, "elec-ex-corr $(scf.xc)", "\n")
    end
end

function make_wannier_centers(scf::self_consistent_field; perturbation=10, norbitals=1)
    centers = []
    ionpos_scf = scf.ionpos
    for ion in ionpos_scf.ionpos
        coords = ion[3:5]
        for i in 1:norbitals
            push!(centers, coords + rand(Float64, 3)/perturbation)
        end
    end
    return centers
end

function write_wannier(wannier::wannier_interpolation, filename::String, scf_filename::String)
    open(filename, create=true, write=true, append=false) do io
        write(io, "include $(scf_filename) \n")
        write(io, "wannier\\ \n")
        write(io, "innerWindow $(wannier.innerWindow[1])   $(wannier.innerWindow[2])\\ \n  ")
        write(io, "outerWindow $(wannier.outerWindow[1])   $(wannier.outerWindow[2])\\ \n  ")
        write(io, "saveWfnsRealSpace", "$(wannier.saveWFNs==true ? "   yes" : "   no")", "\n")
        if wannier.phonon==true
            write(io, "\\ \n")
            for phonon_mesh in wannier.phononSupercell
                write(io, "$(phonon_mesh)  ")
            end
            write(io, "\n")
        end
        write(io, "wannier-initial-state   ", "$(string(scf_filename, ".", "\$", "VAR"))", "\n")
        write(io, "wannier-dump-name   ", "$(string(filename, ".", "\$", "VAR"))", "\n" )
        for center in wannier.wannier_centers
            write(io, "wannier-center Gaussian   ")
            for coord in center
                write(io, "$(coord)", "   ")
            end
            write(io, "\n")
        end
        write(io, "wannier-minimize niterations  $(wannier.wannier_minimize)")
    end
end

function write_phonon(phonon_params::phonon, filename::String, scf_filename::String)
    open(filename, create=true, write=true, append=false) do io
        write(io, "include  $(scf_filename)\n")
        write(io, "initial-state   ", "$(string(scf_filename, ".", "\$", "VAR"))", "\n")
        write(io, "dump-only \n\n")
        write(io, "phonon supercell  ", "$(phonon_params.supercell[1])  $(phonon_params.supercell[2])  $(phonon_params.supercell[3]) ")
    end
end