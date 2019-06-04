using Laplacians
# using DataStructures

# #========================================================================
#      DATA STRUCTURES
# =#


# mutable struct fastQueue
#     q::Vector{Int64}
#     n::Int64
#     curPtr::Int64
#     endPtr::Int64
# end

# fastQueue(n::Int) = fastQueue(zeros(Int64,n), n, 1, 0)

# hasMore(fq::fastQueue) = fq.curPtr <= fq.endPtr

# import Base.push!

# function push!(fq::fastQueue, i)
#     @assert fq.endPtr < fq.n
    
#     fq.endPtr = fq.endPtr + 1
#     fq.q[fq.endPtr] = i
       
# end

# function pull!(fq::fastQueue)
#     @assert hasMore(fq)
    
#     i = fq.q[fq.curPtr]
#     fq.curPtr += 1
    
#     return i
# end

# function reset!(fq::fastQueue)
#     fq.curPtr = 1
#     fq.endPtr = 0
# end

# mutable struct Components
#     c::Vector{SparseMatrixCSC}
#     n::Int64
# end

# Components(n::Int) = Components(zeros(SparseMatrixCSC{Tv,Ti},n), n)

"""Return all disconnected components (their subgraph) as sparseMatrixCSC"""
#"""Returns also the mapping of nodes between subgraphs and original graph"""
function allComp(mat:: SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    cv = components(mat)
    nc = maximum(cv)
    nodes = vecToComps(cv)
    comps = Vector{SparseMatrixCSC{Tv,Ti}}(nc) # (undef nc)
    i = 1
    for c in nodes
        comps[i] = mat[c,c]
        i = i + 1
    end
    if i == nc
        printf("Good!")
    end
    return comps, nodes, nc
end

    # sizes = map(x->length(x), cs)
    # jnk, ind = findmax(sizes)
    # return mat[cs[ind],cs[ind]]
