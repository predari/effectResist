include("structures.jl")
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
function calculateCF(cf:: Array{Float64,1}, n)
    
    for i in 1:n
         cf[i]= n/cf[i]
        #println(cf[i])
    end
    
    return cf
end

function calculateCF(cf, n, sizecf)
    
    for i in 1:sizecf
         cf[i]= n/cf[i]
        #println(cf[i])
    end
    
    return cf
end

function calculateCF(cf :: Float64, n)
    if cf == 0.0
        println("WARNING: effective resistance distance cannot be zero!")
        exit()
    end
    return cf = n/cf
end


##############################################################

function compContractLocalDists(C :: Array{Component,1}, nc :: Int64, path2 :: Array{Int64,2}, dist2 :: Array{Float64,2}, sizes :: Array{Int64,1})
    distcomp = zeros(Float64,nc)
   
    for i in 1:nc
        #println("- C[$i] cf-distances=",C[i].distances)
        # TODO: for (j,p) in enumerate(C[i].path)
        for (j,p) in enumerate(path2[i,:])
            if j == i continue; end
            idx = findin(C[j].link, p)
            #println("In comp=$j node-->$idx sum:", sum(C[j].distances[:,idx+1]))
            distcomp[i] += sizes[j]*dist2[i,j] + sum(C[j].distances[:,idx+1])
        end
        #### following code for debug purposes
        # println(dist2)
        # for (j,p) in enumerate(path2[i,:])
        #     if j == i continue; end
        #     println(sizes[j],"*",dist2[i,j])
        #     distcomp[i] += sizes[j]*dist2[i,j]
        # end
        # println("distcomp after size:", distcomp)
        # for (j,p) in enumerate(path2[i,:])
        #     if j == i continue; end
        #     idx = findin(C[j].link, p)
        #     println("In comp=$j node-->$p sum:", sum(C[j].distances[:,idx+1]))
        #     distcomp[i] += sum(C[j].distances[:,idx+1])
        # end
        # println("distcomp after sum:", distcomp)
    end
    return distcomp
end


function updateLocalDists(C :: Array{Component,1}, sizes :: Array{Int64,1},  edges :: Array{Int64,1}, cmplist :: Array{Int64,1})  
    s :: Int64 = sum(sizes)
    l :: Int64 = length(cmplist)
    for (idx, c) in enumerate(C)
        if size(c.distances,2) < 2
            println("WARNING: # of columns for distances should be at least 2")
        end
        if size(c.distances,2) == 2
            #println("In comp=$idx link with idx=2, (mult with ", s - sizes[idx], ")")
            c.distances[:,2] += (c.distances[:,2] * (s - sizes[idx]))
        else
            for i in 2:size(c.distances,2)
                node = c.link[i-1]
                x = getindex(findin(edges,node))
                if rem.((l-x),2) == 1
                    y = x + 1
                else # (==0)
                    y = x - 1
                end
                yi:: Int64 = cmplist[y]
                #println("In comp=$idx link with idx=$i, (mult with ", sizes[yi], ")")
                c.distances[:,i] += (c.distances[:,i] * sizes[yi])
            end
        end
    end
end


function aggregateLocalDists(C :: Array{Component,1}, distcomp :: Array{Float64,1})
    fdistance = []
    fnodes = []
    for (idx, c) in enumerate(C)
        c.distances = sum(c.distances,2) + distcomp[idx]
        #println("distances:")
        #println(c.distances)
        fdistance = [fdistance ; c.distances]
        fnodes = [fnodes ; c.nodemap]
    end
    return fdistance, fnodes
end

function aggregateDistances(C :: Array{Component,1}, nc :: Int64, path2 :: Array{Int64,2}, dist2 :: Array{Float64,2}, sizes :: Array{Int64,1}, edges :: Array{Int64,1}, cmplist :: Array{Int64,1})
    distcomp = zeros(Float64,nc)
    fdistance = []
    fnodes = []
    for i in 1:nc
        #println("- C[$i] cf-distances=",C[i].distances)
        # TODO: for (j,p) in enumerate(C[i].path)
        for (j,p) in enumerate(path2[i,:])
            if j == i continue; end
            idx = findin(C[j].link, p)
            #println(j," ",p," ",idx)
            distcomp[i] += sizes[j]*dist2[i,j] + sum(C[j].distances[:,idx+1])
        end
    end
    println("distcomp:")
    println(distcomp)
    #println(B.edges," ",length(B.edges)," ",B.comp)
    #sizes -= 1
    for i in 1:2:length(edges)
        #print("(",B.edges[i],",",B.edges[i+1],") \n")
        for (idx, c) in enumerate(cmplist[i:i+1])
            #println(idx," ", c," ", C[c].link)
            for i in 1:length(C[c].link)                   
                #println(B.comp[Int(2/idx)], ": ", sizes[B.comp[Int(2/idx)]])
                C[c].distances[:,i+1] += (C[c].distances[:,i+1] * sizes[ cmplist[Int(2/idx)]])
            end
        end
    end
    for (idx, c) in enumerate(C)
        c.distances = sum(c.distances,2) + distcomp[idx]
        println("distances:")
        println(c.distances)
        fdistance = [fdistance ; c.distances]
        fnodes = [fnodes ; c.nodemap]
    end
    return fdistance, fnodes
end

