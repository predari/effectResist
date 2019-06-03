
function calculateCommuteDists(Linv, n, u) # n is the size of G not L. size(L,1) = n-1
    distances = Array{Array{Float64, 1}}(n-1)
    for i in indices(distances,1) distances[i] = [] end
    for j in 1:n-1
        push!(distances[u],Linv[j,j]) 
    end
    for i in 1:n-1
        for j in i+1:n-1
            push!(distances[i+1], Linv[j,j] + Linv[i,i] - 2*Linv[i,j])
        end
    end
    return distances
end

function calculateNormFullDists(er :: Array{Float64}, n) # n is the size of G not L. size(L,1) = n-1
    distances = Array{Array{Float64, 1}}(n-1)
    for i in indices(distances,1) distances[i] = [] end
    for i in 1:n-1
        for j in i+1:n
            push!(distances[i], (er[i] - er[j])^2)
        end
    end
    return distances
end



function calculateNodeDists(distances, len)
    #n = G.n
    n = len
    cf = zeros(n)
    for i in 1:n-1
        cf[i] = sum(distances[i])
        cf[n] += distances[i][end]
        if i > 1
            for j in i-1:-1:1
                cf[i] += distances[j][i-j]
            end
        end
    end
    return cf
end

#1 ./ [1, 2, 3]
function calculateCF(cf, n)
    
    for i in 1:n
         cf[i]= n/cf[i]
        #println(cf[i])
    end
    
    return cf
end
