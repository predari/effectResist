using Laplacians

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
    tokeep = Array{Ti,1}()
    sizes = zeros(Ti,nc)
    for i in compvec
        sizes[i] += 1
    end
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
    idx = [idx; tokeep]
    l = length(idx)
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
    comps = Vector{SparseMatrixCSC{Tv,Ti}}(nc) # (undef nc)
    for (i, c) in enumerate(nodes)
        comps[i] = mat[c,c]
    end
    return comps, nodes, nc
end
