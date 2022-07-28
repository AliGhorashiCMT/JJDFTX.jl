
function gvectors(filebase::AbstractString)
    iGarr = Vector{Vector{Vector{Float64}}}()
    k = Vector{Vector{Float64}}()
    wk = []
    iGk = Vector{Vector{Float64}}()
    for line in readlines("$filebase.Gvectors")
        #println("gvector")
        parsed_line = filter(x -> !isequal(x, "]"), string.(split(line)))
        parsed_line = filter(x -> !isequal(x, "["), parsed_line)
        if (length(parsed_line) > 0)
            if first(parsed_line) == "#"
                iGk = Vector{Vector{Float64}}()
                push!(k, [parse(Float64, s) for s in parsed_line[3:5]])
                push!(wk, last(parsed_line))
            elseif length(parsed_line) > 1
                push!(iGk, parse.(Float64, parsed_line))
            end
        else
            push!(iGarr, iGk)
        end
    end
    return k, wk, iGarr
end

function return_cg(filebase::AbstractString, kvector::Vector{<:Real}, band::Integer, numbands::Integer)
    ks, _ , iGarr = gvectors(filebase);
    k_idx = findfirst(x -> isapprox(x, kvector, atol=1e-3), ks) 
    println(k_idx)
    wfns = np.fromfile("$filebase.wfns", dtype=np.complex128)
    start_idx = sum(length.(iGarr)[1:k_idx-1])
    num_cgs = length(iGarr[k_idx])
    np.reshape(wfns[1+start_idx*numbands:(start_idx+num_cgs)*numbands], (numbands, num_cgs))[band, :]
end

function return_cg(filebase::AbstractString, ks::Vector{<:Vector{<:Real}},
    iGarr::Vector{<:Vector{<:Vector{<:Real}}}, kvector::Vector{<:Real}, 
    band::Integer, numbands::Integer, spin::Union{Val{'u'}, Val{'d'}})
    
    k_idxs = findall(x -> isapprox(x, kvector, atol=1e-3), ks) 
    k_idx = isa(spin, Val{'u'}) ? first(k_idxs) : k_idxs[2]
    wfns = np.fromfile("$filebase.wfns", dtype=np.complex128)
    start_idx = sum(length.(iGarr)[1:k_idx-1])
    num_cgs = length(iGarr[k_idx])
    np.reshape(wfns[1+start_idx*numbands:(start_idx+num_cgs)*numbands], (numbands, num_cgs))[band, :]
end

function return_cg(filebase::AbstractString, kvector::Vector{<:Real}, band::Integer, numbands::Integer,
    spin::Union{Val{'u'}, Val{'d'}})
    ks, _ , iGarr = gvectors(filebase);
    return_cg(filebase, ks, iGarr, kvector, band, numbands, spin)
end

function real_space_wfn(Ckb::Vector{<:ComplexF64}, iGk::Matrix{<:Real}, S::Vector{<:Integer})
    iGk += np.where(iGk .< 0, reshape(np.repeat(S, size(iGk)[1]), (size(iGk)[1], 3)), 0) 
    stride = [S[3]*S[2], S[3], 1]
    index = Int.(np.dot(iGk, stride)) 
    psi_kb = np.zeros(np.prod(S), dtype=np.complex128)
    psi_kb[index .+ 1] = Ckb
    psi_kb = np.reshape(psi_kb, S) 
    return np.fft.fftn(psi_kb)
end

function return_transformed_gk(gk::Vector{<:Real}, Gvectors::Vector{<:Vector{<:Float64}}, 
    transformation::Array{<:Float64, 2})
    GMatrix = hcat(Gvectors...)
    B = round.(Integer, inv(GMatrix)*transformation*GMatrix)
    return B*gk
end

function kinetic_energy(lattice_vectors::Vector{<:Vector{<:Float64}}, iGk::Vector{<:Vector{<:Float64}}, 
    k::Vector{<:Float64}, Ckb::Vector{<:ComplexF64})

    V = unit_cell_volume(lattice_vectors)/bohrtoangstrom^3
    Gplusksquareds = [sum(a.^2) for a in 
    unnormalize_kvector.(Ref(lattice_vectors), [collect(row) + k for row in eachrow(np.array(iGk))])*bohrtoangstrom];
    sum(((abs.(Ckb)).^2).*Gplusksquareds/2)*V
end