include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("compDistances.jl")
include("appxInvTrace.jl")
include("Linvdiag.jl")

include("bridge.jl")
include("LinvDistance.jl")

using Laplacians
using LightGraphs


function delnodes(A:: SparseMatrixCSC{Float64}, v::Set{Int64})
    idx = setdiff(Array(1:A.n),v)
    return A[idx, idx]
end

# TODO: remove
function delnode2(L, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function erJLT(A:: SparseMatrixCSC{Float64},L:: SparseMatrixCSC{Float64})

    n = A.n
    start_approx_time = time()        
    er = LinvdiagSS(A;JLfac=20)
    u = indmin(er)
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)
    cf = (n > 2000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    end_approx_time = time()
    println("TOTAL APPROX TIME IS: ", end_approx_time - start_approx_time, "(s)")    
  return u, cf;

end

function approx(A:: SparseMatrixCSC{Float64}, L:: SparseMatrixCSC{Float64}, w :: IOStream)
    logw(w,"****** Running (chinese) approx based on JLT ******")
    u, maxcf = erJLT(A,L)
    logw(w,"\t node with argmax{c(", u, ")} = ", maxcf)
    return u
end

function wrapprox(A:: SparseMatrixCSC{Float64}, L:: SparseMatrixCSC{Float64}, w :: IOStream, solution = nothing)
    u = approx(A, L, w)
    if solution != nothing
        if solution != u
            logw(w,"\t THIS APPROX RESULT IS DIFFERENT THAN OTHERS (OR THE EXACT SOL = ", solution,")")
        end        
    end
end


# TODO: remove
function erINV(G, alldistances)    
    n = G.n;
    er = zeros(n)
    L = lapl(G)
    u = 1
    L2 = delnode2(L,u,n)
    Linv = inv(L2)
    distances = calculateCommuteDists(Linv, n, u)
    return distances
end

function exact(G, w :: IOStream)
    logw(w,"****** Running (my) exact ******")
    distances = erINV(G,1)
    #println("Exact distance:",distances)
    cf = calculateNodeDists(distances, G.n)
    #println("Exact distance:",cf)
    cf = calculateCF(cf, G.n)
    logw(w,"\t node with argmax{c(", indmax(cf), ")} = ",
         maximum(cf))
    return indmax(cf)
end


function approxcore2(A:: SparseMatrixCSC{Float64},L:: SparseMatrixCSC{Float64}, w :: IOStream)
    logw(w,"****** Running (core2) approx ******")
    g = LightGraphs.Graph(A)
    n = A.n
    println("A.n = ",n)
    core2bdry = Array{Int64,1}
    sizes = Array{Array{Int, 1}}
    core2bdry, sizes = bridges3(g)
    core2 = k_core(g, 2)
    #println("size of core2 = ",length(core2))
    if isempty(setdiff(core2bdry,core2)) == false
        println("WARNING: boundary nodes of core2 are not in core2!")
    end    
    distances = LinvDistance(A[core2,core2], core2bdry, sizes, core2)
    cf = calculateCF(distances, n, length(distances))
    logw(w,"\t node with argmax{c(",  core2[indmax(cf)], ")} = ", maximum(cf))
    return core2[indmax(cf)]
end

function wrapproxcore2(A:: SparseMatrixCSC{Float64}, L:: SparseMatrixCSC{Float64}, w :: IOStream, solution = nothing)
    u = approxcore2(A, L, w)
    if solution != nothing
        if solution != u
            logw(w,"\t THIS APPROX RESULT IS DIFFERENT THAN OTHERS (OR THE EXACT SOL = ", solution,")")
        end        
    end
end
