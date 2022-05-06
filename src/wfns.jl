
function gvectors(filebase::AbstractString)
    iGarr = Vector{Vector{Vector{Float64}}}()
    k = Vector{Vector{Float64}}()
    wk = []
    iGk = Vector{Vector{Float64}}()
    for line in readlines("$filebase.Gvectors")
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
    #println(numbands)
    wfns = np.fromfile("$filebase.wfns", dtype=np.complex128)
    start_idx = sum(length.(iGarr)[1:k_idx-1])
    num_cgs = length(iGarr[k_idx])
    np.reshape(wfns[1+start_idx*numbands:(start_idx+num_cgs)*numbands], (numbands, num_cgs))[band, :]
end

function return_transformed_gk(gk::Vector{<:Real}, Gvectors::Vector{<:Vector{<:Float64}}, 
    transformation::Array{<:Float64, 2})
    GMatrix = hcat(Gvectors...)
    B = round.(Integer, inv(GMatrix)*transformation*GMatrix)
    #println(B)
    return B*gk

end