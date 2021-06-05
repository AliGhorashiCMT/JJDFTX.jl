function parsetonum(stringtoparse::AbstractString, typetoparse::Union{Type{Int}, Type{Float64}})
    splittedstring = string.(split(stringtoparse))
    num = 0 
    for split in splittedstring
        try 
            num = parse(typetoparse, split)
        catch
        end
    end
    return num
end

function parsetonums(stringtoparse::AbstractString, typetoparse::Union{Type{Int}, Type{Float64}})
    splittedstring = string.(split(stringtoparse))
    nums = Vector{typetoparse}() 
    for split in splittedstring
        try 
            push!(nums, parse(typetoparse, split))
        catch
        end
    end
    return nums
end


"""
$(TYPEDSIGNATURES)

Plot a 2d slice of the electronic density in a material given an axis perpendicular to the 2d slice.
Note that to obtain a good plot of the density, a high density cutoff (which will make a finer density grid) should be used.
"""
function plot_density(density_file::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, repeatnums::Vector{<:Integer} = [0, 0], slice::Integer=1; kwargs...)
    n = np.fromfile(density_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        contains(r, "Chosen fftbox") || continue
        S = parsetonums(r, Int)
        break ##Only look at first instance of fftbox 
    end
    V=0
    for r in readlines(outfile)
        contains(r, "volume") || continue
        V = parsetonum(r, Float64)
        break #In cases of coulomb truncation, the volume is doubled afterwards, so we only look at the first instance of volume in the output file
    end
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    nelectrons = dV * np.sum(n)
    println("Nelectrons = ", nelectrons )
    n = np.reshape(n, S)
    if perpaxis isa Val{'x'}
        n = n[slice,:,:]        
        n = np.roll(n, Int(S[2]/2), axis=0) 
        n = np.roll(n, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        n = n[:,slice,:]        
        n = np.roll(n, Int(S[1]/2), axis=0) 
        n = np.roll(n, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        n = n[:,:,slice]        
        n = np.roll(n, Int(S[1]/2), axis=0) 
        n = np.roll(n, Int(S[2]/2), axis=1) 
    end
    nprime = n
    for _ in 1:repeatnums[1]
        nprime = hcat(nprime, n)
    end
    n = nprime
    for _ in 1:repeatnums[2]
        nprime = vcat(nprime, n)
    end
    display(heatmap(nprime; kwargs...))
    return V, S, nelectrons
end

plot_density(density_file::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, repeatnum::Integer=0, slice::Integer=1; kwargs...) = plot_density(density_file, outfile, perpaxis, [repeatnum, repeatnum], slice; kwargs...)

"""
$(TYPEDSIGNATURES)

Plot a 2d slice of the spin density in a material. This assumes that a spin polarized calculation has been done for which spin polarized densities have been dumped. 
This is particularly useful in situations where one is interested, for instance, in the hole localization of a polaron or the spatial properties of a defect state. 
"""
function plot_diffdensity(density_file1::AbstractString, density_file2::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, repeatnums::Vector{<:Integer}=[0, 0], slice::Integer=1; kwargs...)
    nup = np.fromfile(density_file1, dtype=np.float64)   
    ndn = np.fromfile(density_file2, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        contains(r, "Chosen fftbox") || continue
        S = parsetonums(r, Int)
        break ##Only look at first instance of fftbox 
    end
    V=0
    for r in readlines(outfile)
        contains(r, "volume") || continue
        V = parsetonum(r, Float64)
        break #In cases of coulomb truncation, the volume is doubled afterwards, so we only look at the first instance of volume in the output file
    end
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    numelectronsup = dV * np.sum(nup)
    numelectronsdn = dV * np.sum(ndn)

    println("Nelectrons of Species 1= ", numelectronsup)
    println("Nelectrons of Species 2= ", numelectronsdn)

    nup = np.reshape(nup, S)
    ndn = np.reshape(ndn, S)

    if perpaxis isa Val{'x'}
        nup = nup[slice,:,:]        
        nup = np.roll(nup, Int(S[2]/2), axis=0) 
        nup = np.roll(nup, Int(S[3]/2), axis=1) 
        ndn = ndn[slice,:,:]        
        ndn = np.roll(ndn, Int(S[2]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        nup = nup[:,slice,:]        
        nup = np.roll(nup, Int(S[1]/2), axis=0) 
        nup = np.roll(nup, Int(S[3]/2), axis=1)
        ndn = ndn[:,slice,:]        
        ndn = np.roll(ndn, Int(S[1]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        nup = nup[:,:,slice]        
        nup = np.roll(nup, Int(S[1]/2), axis=0) 
        nup = np.roll(nup, Int(S[2]/2), axis=1) 
        ndn = ndn[:,:,slice]        
        ndn = np.roll(ndn, Int(S[1]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[2]/2), axis=1) 
    end

    ndiff = nup-ndn
    nprime = ndiff
    for _ in 1:repeatnums[1]
        nprime = hcat(nprime, ndiff)
    end
    ndiff = nprime
    for _ in 1:repeatnums[2]
        nprime = vcat(nprime, ndiff) 
    end

    display(heatmap(nprime; kwargs...))
    return V, S, numelectronsup, numelectronsdn
end

plot_diffdensity(density_file1::AbstractString, density_file2::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, repeatnum::Integer=0, slice::Integer=1; kwargs...) = plot_diffdensity(density_file1, density_file2, outfile, perpaxis, [repeatnum, repeatnum], slice; kwargs...)

##Plots the wavefunction density 
function plot_wfns(wfn_file::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, component::Union{Val{'r'}, Val{'i'}, Val{'a'}} )
    @warn "For accurate results, make sure norm conserving pseudopotentials are being used"
    ψ = np.fromfile(wfn_file, dtype=np.complex)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        contains(r, "Chosen fftbox") || continue
        if length(S) != 0 ##The wavefunction grid always comes second in the output file, so reset the grid if already read
            S = Vector{Int}()
        end
        splittedfft = split(r)
        for _ in splittedfft
            try
                a = parse(Int, split)
                push!(S, a)
            catch

            end
        end
    end
    V=0
    for r in readlines(outfile)
        contains(r, "volume") || continue
        splittedV = split(r)
        for _ in splittedV
            try 
                V = parse(Float64, split)
            catch

            end
        end
    end
    println("The wavefunction grid is: ", S)
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    ψsquared = np.abs(ψ).^2
    ψimag = imag(ψ)
    ψreal = real(ψ)
    println("Normalization of wavefunction: ", dV * np.sum(ψsquared))
    ψsquared = np.reshape(ψsquared, S)
    ψimag = np.reshape(ψimag, S)
    ψreal = np.reshape(ψreal, S)
    if perpaxis isa Val{'x'}
        ψsquared = ψsquared[1,:,:]        
        ψsquared = np.roll(ψsquared, Int(S[2]/2), axis=0) 
        ψsquared = np.roll(ψsquared, Int(S[3]/2), axis=1) 
        ψreal = ψreal[1,:,:]        
        ψreal = np.roll(ψreal, Int(S[2]/2), axis=0) 
        ψreal = np.roll(ψreal, Int(S[3]/2), axis=1) 
        ψimag = ψimag[1,:,:]        
        ψimag= np.roll(ψimag, Int(S[2]/2), axis=0) 
        ψimag = np.roll(ψimag, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        ψsquared = ψsquared[:,1,:]        
        ψsquared = np.roll(ψsquared, Int(S[1]/2), axis=0) 
        ψsquared = np.roll(ψsquared, Int(S[3]/2), axis=1) 
        ψreal = ψreal[:,1,:]        
        ψreal = np.roll(ψreal, Int(S[1]/2), axis=0) 
        ψreal = np.roll(ψreal, Int(S[3]/2), axis=1) 
        ψimag = ψimag[:,1,:]        
        ψimag = np.roll(ψimag, Int(S[1]/2), axis=0) 
        ψimag = np.roll(ψimag, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        ψsquared = ψsquared[:,:,1]        
        ψsquared = np.roll(ψsquared, Int(S[1]/2), axis=0) 
        ψsquared = np.roll(ψsquared, Int(S[2]/2), axis=1) 
        ψreal = ψreal[:,:,1]        
        ψreal= np.roll(ψreal, Int(S[1]/2), axis=0) 
        ψreal = np.roll(ψreal, Int(S[2]/2), axis=1) 
        ψimag = ψimag[:,:,1]        
        ψimag = np.roll(ψimag, Int(S[1]/2), axis=0) 
        ψimag = np.roll(ψimag, Int(S[2]/2), axis=1) 
    end
    if component isa Val{'i'}
        heatmap(ψimag)
    elseif component isa Val{'r'}
        heatmap(ψreal)
    elseif component isa Val{'a'}
        heatmap(ψsquared)
    end
end

"""
$(TYPEDSIGNATURES)

Compuate the wavefunction overlap from dumped wavefunction files. This function is provided mostly as a check of normalization. 

"""
function wavefunctionoverlap(wfn_file1::AbstractString, wfn_file2::AbstractString, outfile::AbstractString)
    @warn "Make sure norm conserving (e.g. SG15 Pseudopotentials) are being used\n\n"
    ψ₁ = np.fromfile(wfn_file1, dtype=np.complex) 
    ψ₂ = np.fromfile(wfn_file2, dtype=np.complex)
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        contains(r, "Chosen fftbox") || continue
        if length(S) != 0 ##The wavefunction grid always comes second in the output file, so reset the grid if already read
            S = Vector{Int}()
        end
        splittedfft = split(r)
        for _ in splittedfft
            try
                a = parse(Int, split)
                push!(S, a)
            catch
            end
        end
    end
    V=0
    for r in readlines(outfile)
        contains(r, "volume") || continue
        splittedV = split(r)
        for _ in splittedV
            try 
                V = parse(Float64, split)
            catch
            end
        end
    end
    println("The wavefunction grid is: ", S)
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    overlap = sum(conj(ψ₁).*ψ₂)*dV
    println("The overlap of the two provided wavefunctions is: ", overlap)
    #return overlap
end

"""
$(TYPEDSIGNATURES)

Plot a 2d slice of the self consistent potential. 
"""
function plot_scfpotential(scf_file::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, )
    scf = np.fromfile(scf_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        contains(r, "Chosen fftbox") || continue
        splittedfft = split(r)
        for _ in splittedfft
            try
                a = parse(Int, split)
                push!(S, a)
            catch

            end
        end
        break ##Only look at first instance of fftbox 
    end
    scf = np.reshape(scf, S)
    if perpaxis isa Val{'x'}
        scf = scf[1,:,:]        
        scf = np.roll(scf, Int(S[2]/2), axis=0) 
        scf = np.roll(scf, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        scf = scf[:,1,:]        
        scf = np.roll(scf, Int(S[1]/2), axis=0) 
        scf = np.roll(scf, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        scf = scf[:,:,1]        
        scf = np.roll(scf, Int(S[1]/2), axis=0) 
        scf = np.roll(scf, Int(S[2]/2), axis=1) 
    end
    heatmap(scf)
end

"""
$(TYPEDSIGNATURES)

Plot a 2d slice of the electrostatic potential. 
"""
function plot_dtot(dtot_file::AbstractString, outfile::AbstractString, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}} )
    dtot = np.fromfile(dtot_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        contains(r, "Chosen fftbox") || continue
        splittedfft = split(r)
        for _ in splittedfft
            try
                a = parse(Int, split)
                push!(S, a)
            catch

            end
        end
        break ##Only look at first instance of fftbox 
    end
    dtot = np.reshape(dtot, S)
    if perpaxis isa Val{'x'}
        dtot = dtot[1,:,:]        
        dtot = np.roll(dtot, Int(S[2]/2), axis=0) 
        dtot = np.roll(dtot, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        dtot = dtot[:,1,:]        
        dtot = np.roll(dtot, Int(S[1]/2), axis=0) 
        dtot = np.roll(dtot, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        dtot = dtot[:,:,1]        
        dtot = np.roll(dtot, Int(S[1]/2), axis=0) 
        dtot = np.roll(dtot, Int(S[2]/2), axis=1) 
    end
    heatmap(dtot)
end