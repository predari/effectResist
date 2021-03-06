using Laplacians

# function components(mat::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
#   n = mat.n

#   order = Array{Ti}(undef, n)
#   comp = zeros(Ti,n)

#   # note that all of this casting is unnecessary.
#   # but, some of it speeds up the code
#   # I have not figured out the minimal necessary
#   c::Ti = 0

#   colptr = mat.colptr
#   rowval = mat.rowval

#   @inbounds for x in 1:n
#     if (comp[x] == 0)
#       c = c + 1
#       comp[x] = c

#       if colptr[x+1] > colptr[x]
#         ptr::Ti = 1
#         orderLen::Ti = 2
#         order[ptr] = x

#         while ptr < orderLen
#           curNode = order[ptr]

#           for ind in colptr[curNode]:(colptr[curNode+1]-1)
#             nbr = rowval[ind]
#             if comp[nbr] == 0
#               comp[nbr] = c
#               order[orderLen] = nbr
#               orderLen += 1
#             end # if
#           end # for
#           ptr += 1
#         end # while
#       end # if
#     end

#   end

#   return comp
# end # function


function nodesInComps(compvec::Vector{Ti}) where Ti
    nc = maximum(compvec)
    sizes = zeros(Ti,nc)
    for i in compvec
        sizes[i] += 1
    end
    idx = findin(sizes,1)
    # do not create components for comps of size 1
    #deleteat!(compvec,idx)
    # for i in idx
    #     deleteat!(compvec,findfirst(compvec,i))
    # end
    # #filter!(x->x \nin 10,a)
    # deleteat!(sizes,idx)

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
    ### TODO: do that before constructing comps Vector!!
    deleteat!(comps,idx)
    #println(comps)

    return comps
end # vecToComps


# Components(n::Int) = Components(zeros(SparseMatrixCSC{Tv,Ti},n), n)

"""Return all disconnected components (their subgraph) as sparseMatrixCSC"""
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

    # sizes = map(x->length(x), cs)
    # jnk, ind = findmax(sizes)
    # return mat[cs[ind],cs[ind]]
