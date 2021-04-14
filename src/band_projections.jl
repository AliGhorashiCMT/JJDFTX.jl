function bandprojections(filebase::String, nband::Integer, norbital::Integer, kindex::Integer, nbands::Integer)
    np.loadtxt("$(filebase).bandProjections", skiprows=3)[(kindex-1)*nbands+nband, norbital]
end

function plotbandprojections(filebase::String, nband::Integer, nbands::Integer)
    projections = []
    numorbitals = size(np.loadtxt("$(filebase).bandProjections", skiprows=3))[2]
    for o in 1:numorbitals
        projectionperk = []
        for k in 1:100
            try
                push!(projectionperk, bandprojections(filebase, nband, o, k, nbands))
            catch
                break
            end
        end
        push!(projections, projectionperk)
    end
    plot()
    plot!(projections)
end