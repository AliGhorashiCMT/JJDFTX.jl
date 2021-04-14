function plot_density(density_file::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}})
    n = np.fromfile(density_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
            break ##Only look at first instance of fftbox 
        end
    end
    V=0
    for r in readlines(outfile)
        if contains(r, "volume")
            splittedV = split(r)
            for split in splittedV
                try 
                    V = parse(Float64, split)
                catch

                end
            end
        end
    end
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    println("Nelectrons = ", dV * np.sum(n))
    n = np.reshape(n, S)
    if perpaxis isa Val{'x'}
        n = n[1,:,:]        
        n = np.roll(n, Int(S[2]/2), axis=0) 
        n = np.roll(n, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        n = n[:,1,:]        
        n = np.roll(n, Int(S[1]/2), axis=0) 
        n = np.roll(n, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        n = n[:,:,1]        
        n = np.roll(n, Int(S[1]/2), axis=0) 
        n = np.roll(n, Int(S[2]/2), axis=1) 
    end
    heatmap(n)
end

function plot_diffdensity(density_file1::String, density_file2::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}})
    nup = np.fromfile(density_file1, dtype=np.float64)   
    ndn = np.fromfile(density_file2, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
            break ##Only look at first instance of fftbox 
        end
    end
    V=0
    for r in readlines(outfile)
        if contains(r, "volume")
            splittedV = split(r)
            for split in splittedV
                try 
                    V = parse(Float64, split)
                catch

                end
            end
        end
    end
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    println("Nelectrons of Species 1= ", dV * np.sum(nup))
    println("Nelectrons of Species 2= ", dV * np.sum(ndn))
    nup = np.reshape(nup, S)
    ndn = np.reshape(ndn, S)
    if perpaxis isa Val{'x'}
        nup = nup[1,:,:]        
        nup = np.roll(nup, Int(S[2]/2), axis=0) 
        nup = np.roll(nup, Int(S[3]/2), axis=1) 
        ndn = ndn[1,:,:]        
        ndn = np.roll(ndn, Int(S[2]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        nup = nup[:,1,:]        
        nup = np.roll(nup, Int(S[1]/2), axis=0) 
        nup = np.roll(nup, Int(S[3]/2), axis=1)
        ndn = ndn[:,1,:]        
        ndn = np.roll(ndn, Int(S[1]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        nup = nup[:,:,1]        
        nup = np.roll(nup, Int(S[1]/2), axis=0) 
        nup = np.roll(nup, Int(S[2]/2), axis=1) 
        ndn = ndn[:,:,1]        
        ndn = np.roll(ndn, Int(S[1]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[2]/2), axis=1) 
    end
    heatmap(nup-ndn)
end

##Plots the wavefunction density 
function plot_wfns(wfn_file::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, component::Union{Val{'r'}, Val{'i'}, Val{'a'}} )
    @warn "For accurate results, make sure norm conserving pseudopotentials are being used"
    ψ = np.fromfile(wfn_file, dtype=np.complex)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            if length(S) != 0 ##The wavefunction grid always comes second in the output file, so reset the grid if already read
                S = Vector{Int}()
            end
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
        end
    end
    V=0
    for r in readlines(outfile)
    if contains(r, "volume")
        splittedV = split(r)
        for split in splittedV
            try 
                V = parse(Float64, split)
            catch

            end
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

function wavefunctionoverlap(wfn_file1::String, wfn_file2::String, outfile::String)
    @warn "Make sure norm conserving (e.g. SG15 Pseudopotentials) are being used\n\n"
    ψ₁ = np.fromfile(wfn_file1, dtype=np.complex) 
    ψ₂ = np.fromfile(wfn_file2, dtype=np.complex)
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
    if contains(r, "Chosen fftbox")
        if length(S) != 0 ##The wavefunction grid always comes second in the output file, so reset the grid if already read
            S = Vector{Int}()
        end
        splittedfft = split(r)
        for split in splittedfft
            try
                a = parse(Int, split)
                push!(S, a)
            catch

            end
        end
    end
    end
    V=0
    for r in readlines(outfile)
    if contains(r, "volume")
        splittedV = split(r)
        for split in splittedV
            try 
                V = parse(Float64, split)
            catch

            end
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

function plot_scfpotential(scf_file::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, )
    scf = np.fromfile(scf_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
            break ##Only look at first instance of fftbox 
        end
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

function plot_dtot(dtot_file::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}}, )
    dtot = np.fromfile(dtot_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
            break ##Only look at first instance of fftbox 
        end
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