include("structures.jl")
using Laplacians
using DataStructures ## for SortedSets



## return the numbering of nodes for each component
function nodesInComps(compvec::Vector{Ti}) where Ti

    nc = maximum(compvec)
    sizes = zeros(Ti,nc)
    for i in compvec
        sizes[i] += 1
    end
    idx = findin(sizes,1)
    comps = Vector{Vector{Ti}}(nc)
     for i in 1:nc
         comps[i] = zeros(Ti,sizes[i])
     end

    ptrs = zeros(Ti,nc)
     for i in eachindex(compvec)
        c = compvec[i]
        ptrs[c] += 1
        comps[c][ptrs[c]] = i
     end
    deleteat!(comps,idx)
    return comps
end 

#### do not delete all components of size 1
function nodesInComps(compvec::Vector{Ti}, bridges :: Array{Ti,1}) where Ti
    nc = maximum(compvec)
    #println("compvec", compvec)
    tokeep = Array{Ti,1}()
    sizes = zeros(Ti,nc)
    for i in compvec
        sizes[i] += 1
    end
    #println(sizes)
    ### TODO: maybe count func is slow?? Check
    for u in unique(bridges)
        if count(x->x==u,bridges) > 1
            push!(tokeep,compvec[u])
        end
    end

    ### previous implementation:
    # idx = findin(sizes,1)
    # for u in unique(bridges)
    #     if count(x->x==u,bridges) > 1
    #         push!(tokeep,u)
    #     end
    # end
    # idx2 = setdiff(idx,tokeep)
    idx = Array{Ti,1}()
    for (i,size) in enumerate(sizes)
        if size != 1
            push!(idx,i)
        end
    end
    l = length(idx)
    println("l = $l and idx= ", idx)
    idx = [idx; tokeep]
    l = length(idx)
    println("l = $l and idx= ", idx)
    #### TODO : I just changed that!! Not sure if idx should unique
    idx = unique(idx)
    l = length(idx)
    println("l = $l and idx= ", idx)
    comps = Vector{Vector{Ti}}(l)
     for i in 1:l
         comps[i] = zeros(Ti,sizes[idx[i]])
     end    
    
    ptrs = zeros(Ti,l)
     for i in eachindex(compvec)
         c = compvec[i]
         if in(c,idx) == false
             continue;
         end
         j = getindex(findin(idx,c))
         ptrs[j] += 1
         comps[j][ptrs[j]] = i
     end    
    # println("comps size= ",size(comps,1))
    # for i in 1:size(comps,1)
    #     println("$i ", comps[i])
    # end
    ### previous implementation:
    # deleteat!(comps,idx)
    return comps
end 

    
# Components(n::Int) = Components(zeros(SparseMatrixCSC{Tv,Ti},n), n)
"""Return all disconnected components (their subgraph) as sparseMatrixCSC. components of size one 
are all removed!"""
#"""Returns also the mapping of nodes between subgraphs and original graph"""
function allComp(mat:: SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    cv = components(mat)
    nodes = nodesInComps(cv)
    nc = length(nodes)
    comps = Vector{SparseMatrixCSC{Tv,Ti}}(nc) # (undef nc)
    for (i, c) in enumerate(nodes)
        comps[i] = mat[c,c]
    end
    return comps, nodes, nc
end

"""Return all disconnected components (their subgraph) as sparseMatrixCSC. components of size one 
that belong to bridges are kept!!"""
#"""Returns also the mapping of nodes between subgraphs and original graph"""

function allComp(mat:: SparseMatrixCSC{Tv,Ti}, bridges :: Array{Ti,1}) where {Tv,Ti}
    cv = components(mat)
    nodes = nodesInComps(cv,bridges)
    nc = length(nodes)
    #println("nc:", nc)
    comps = Vector{SparseMatrixCSC{Tv,Ti}}(nc) # (undef nc)
    for (i, c) in enumerate(nodes)
        #println(c)
        comps[i] = mat[c,c]
    end
    return comps, nodes, nc
end


function buildComponents(A :: SparseMatrixCSC{Float64}, B :: Bridges)
    start_time = time()
    cmps, map, ncmps = allComp(A, B.edges)
    println("finding components time: ", time() - start_time, "(s)")
    t = time()
    C = Array{Component,1}(ncmps)
    maxc = 0;
    #### Is : for i = eachindex(a) faster than for i = 1:n?
    for i in eachindex(cmps) 
        bdry = collect(intersect(DataStructures.SortedSet(B.core2nodes), Set(map[i])))
        # following command is faster but doesn't preserve order
        # bdry = collect(intersect(Set(B.core2nodes), Set(map[i])))
        # TODO: check if intersect is still that bad with (sorted) arrays (compared to Sets)
        # link array doesn't need to be sorted.
        link = collect(intersect(Set(B.edges), Set(map[i])))
        index = findin(B.core2nodes,bdry)
        #println(B.ext[index])
        tcount :: Int64 = 0
        for j in B.ext[index]
            for k in j
                tcount += k
            end
        end
        B.comp[findin(B.edges,link)] = i 
        C[i] = Component(cmps[i],cmps[i].n,cmps[i].n + tcount, map[i],bdry,link,length(link),
                         zeros(cmps[i].n,length(link)), B.ext[index])
        if cmps[i].n >= maxc
            maxc = i
        end
    end
    println("building structure components time: ", time() - t, "(s)")
    return C
end

function compRealSizes(C :: Array{Component,1}, nc :: Int64)
    sizes = zeros(Int64,nc)
    for i in 1:nc
        externalNodes = 0
        for j in C[i].external
            for k in j
                externalNodes += k
            end
        end
        sizes[i] = C[i].nc + externalNodes
    end
    return sizes
end


#########   TODO: create checkDistances(dist1 :: Array{Float64,1},dist2 :: Array{Float64,1})
########    to check at which extend distances are similar.


function stripcore1nodes(comp  :: Array{Int64,1} , core3nodes :: Array{Int64,1})
    if length(comp) != length(core3nodes)
        println("WARNING: length of node-vector and their component index vector should be the same")
    end
    pst = []
    for i in eachindex(comp)
        if comp[i] == 0
            push!(pst,i)
        end
    end
    deleteat!(comp,pst)
    deleteat!(core3nodes,pst)
    return comp, core3nodes
end


function contractAdjGraph(edges :: Array{Int64,1}, cmplist :: Array{Int64,1} , C :: Array{Component,1}, nc :: Int64)
    
    ue = unique(edges)
    n = length(ue)
    #println(edges)
    #println(ue)
    #println(1:n)
    #println(cmplist)
    ## mape mapc store the component index and the node index
    ## for the numbering of the contracted graph.
    mape = zeros(Int64,n)
    mapc = zeros(Int64,n)
    cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
    for i in 1:2:length(edges)   
        x = getindex(findin(ue,edges[i]))
        y = getindex(findin(ue,edges[i+1]))
        mapc[x] = cmplist[i]
        mapc[y] = cmplist[i+1]
        mape[x] = edges[i]
        mape[y] = edges[i+1]
        #cA[i,i+1] .= 1.0
        #cA[i+1,i] .= 1.0
        cA[x,y] .= 1.0
        cA[y,x] .= 1.0

    end
    #println(mapc)
    #println(mape)
    #println(full(cA))
    ### TODO:: address the fact that there are components with size of 1.
    ### if this is not addressed, then we have edges that are not bridges
    ### for instance 1--C0---25, 1--C'0---4, which after the strip function
    ### will result to 1--25 and 1--4, which is not allowed and creates
    ### non unique edge indexes that result in different of length
    ### between idx and lidx2, lidx1 !
    for i in 1: nc
        idx = findin(mapc,i)
        #println("i $i ", idx)
        if length(idx) > 1
            c = C[i]
            ### TODO: I don't really need the following if
            if c.nc == 1
                println("WARNING: component size should be larger than one!")
            end
            #println("component:",i)
            #println(idx," ", edges[idx])
            #println("c.link size = ", length(c.link))
            #println(edges[idx])
            #println(c.link)
            #println(setdiff(c.link,edges[idx]))
            #lidx1 = findin(c.link, edges[idx])
            lidx1 = findin(c.link, mape[idx])
            #println(" ", length(lidx1))
            #lidx2 = findin(c.nodemap, edges[idx])
            lidx2 = findin(c.nodemap, mape[idx])
            lidx1 += 1
            #println(lidx1, lidx2)
            dist = c.distances[lidx2,lidx1]
            #println(lidx2)
            #println(lidx1)
            #println("idx = ",  " len= ", length(idx), " ", length(idx))
            #println("dist = ", " len= ", size(dist,1), " ",size(dist,2))
            #len = length(idx)
            cA[idx,idx] = dist
        end
    end
    ### need to clear ones in diagonal caused by cA[idx,idx] = dist
    for i in 1:n
        if cA[i,i] != 0.0
            cA[i,i] = 0.0
        end
    end
    dropzeros!(cA)
    
    if !Laplacians.isConnected(cA)
        logw(w,"WARNING: contracted graph should be connected! ");
        logw(w,"WARNING: Program will exit!");
        exit()
    end
    # if !(isConnected(cA) && nnz(cA) == 2*(cA.n-1))
    #     logw(w,"WARNING: contracted graph should be a tree! ");
    # end
    return cA
end

