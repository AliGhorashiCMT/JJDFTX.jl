
#The methods here are designed after "Plasmonics in Argentene by Shankar"

#Below, we try to calculate 1/τ(ω) as defined in cited paper above. 

function tauinverse(e::Real, μ::Real; histogram_width::Integer=10, mesh::Integer=100)
    tauinv = 0 
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ek = 6*k 
        abs(ek-μ)*histogram_width < 0.5 || continue
        for (xmesh1, ymesh1) in Tuple.(CartesianIndices(rand(mesh, mesh)))

        end
        tauinv += histogram_width/(mesh^2)

    end
    return (2π/(ħ*gμ))*summand
end