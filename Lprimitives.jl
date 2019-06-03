include("graph.jl")
include("lapl.jl")
include("logw.jl")
#include("bridge.jl")
#include("components.jl")
#include("Exact.jl")
include("appxInvTrace.jl")
#include("Linvdiag.jl")
include("compDistances.jl")
using Laplacians
using LightGraphs

# function erJLT(a::SparseMatrixCSC{Tv,Ti}, edgevalues; ep=0.3, matrixConcConst=4.0, JLfac=200.0) where {Tv,Ti}
#     n = size(a,1)
#     # time start
#     start_time = time()
#     f = approxCholLap(a,tol=1e-5);
    
#     println(f)
#     k = round(Int, JLfac*log(n)) # number of dims for JL
#     # U = incident matrix: UU* = L
#     U = wtedEdgeVertexMat(a)
#     m = size(U,1)
#     println("m = ",m, ", n = ", n)
#     cf = zeros(n)
#     er_edge = zeros(n)
#     er_edge_table = zeros(n,n)
#     println("q projections = ",k)
#     k = 2
#     for i = 1:k
#         r = randn(m)
#         println("vector r")
#         println(r)
#         ur = U'*r
#         println("vector ur")
#         println(ur)
#         v = zeros(n)
#         v = f(ur[:])
#         println("vector v")
#         println(v)
#         cf.+= v.^2/k
#         for i = 1:n
#             for j = 1:n    
#                 er_edge_table[i,j] = er_edge_table[i,j] + norm(v[i]-v[j],2)/k
#             end 
#         end
#         println("vector er")
#         println(cf)
#     end
#     println(er_edge_table)
#     cf_edge = sum(er_edge_table, 2)
#     println("cf_edge:")
#     println(cf_edge)
#     println("cf:")
#     println(cf)
#     # time finish
#     elapsed_time = time()-start_time
#     print(" * erJLT total time: ",time() - elapsed_time, " (s)\n")
    
#   return cf;
# end


# see LinvdiagSS
function LinvSS(a::SparseMatrixCSC{Float64}; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    f = approxCholLap(a,tol=1e-5);

    n = size(a,1)
    println("n=",n)
    k = round(Int, JLfac*log(n)) # number of dims for JL

    U = wtedEdgeVertexMat(a)
    m = size(U,1)
    println("m=",m)
    R = randn(m,k)
    UR = U'*R;

    V = zeros(n,k)
    k = 2
    for i = 1:k
        V[:,i] .= f(UR[:,i])
        #println(typeof(V))
        #println(size(V,1),",",size(V,2))        
        #println(typeof(r))
        #V[:,i] .= V[:,i]*R[:,i]
        #rst += sum(v.*r)
    end
    return V/k
end

function LinvdiagEdgeDist(a::SparseMatrixCSC{Float64}; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    f = approxCholLap(a,tol=1e-5);
    
    n = size(a,1)
    k = round(Int, JLfac*log(n)) # number of dims for JL
    
    U = wtedEdgeVertexMat(a)
    m = size(U,1)
    er = zeros(n)
    dr = zeros(n)
    cf = zeros(n)
    
    distances = Array{Array{Float64, 1}}(n-1)
    for i in indices(distances,1) distances[i] = [] end
    # for i in 1:n-1
    #     for j in i+1:n
    #         push!(distances[i], 0.0)
    #     end
    # end
    for i = 1:k # q 
        # random gaussian projections Q
        r = randn(m) 
        # compute (QW^(1/2)B)
        ur = U'*r 
        v = zeros(n)
        # solve(L, ((QW^(1/2)B)^T))^T
        # v here is Z_{i,:}
        v = f(ur[:])
        #println("vector v=",v)
        #println("type=",typeof(v),"size=",size(v,1),",",size(v,2))
        # er(u) = Sigma_{v\in V\u}(Z_{u,v})^2
        # decause I'll do the above k times
        # given the k random projections,
        # I have to divide by k
        dr.+= v./k
        er.+= v.^2/k
    end
    # for i in 1:n-1
    #     for j in 1:n-i
    #         #println(i,",",j,",",j+i)
    #         distances[i][j] = er[i] + er[j+i] - 2*dr[i]*dr[j+1]
    #     end
    # end
    for i in 1:n
        cf[i] = sum(er) + 4*er[i] -2(sum(dr))
    end

    #return distances;
    return cf
end

function erJLT(G)
    n = G.n
    
    A, L = sparseAdja(G)
    er = LinvdiagSS(A;JLfac=20)
    u = indmin(er)
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)
    cf = (n > 200) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
  return u, cf;

end

function erJLT(G, alldistances)
    n = G.n
    cf = zeros(n)
    A, L = sparseAdja(G)
    #er = LinvdiagSS(A;JLfac=20)
    #er = LinvdiagEdgeDist(A;JLfac=20)
    #println(er)
    #u = indmin(er)
    
    #u = 1
    #L2 = delnode2(L,u,n)
    #A2 = delnode2(A,u,n)
    distances = LinvdiagEdgeDist(A;JLfac=20)
    println("distances size:", size(distances,1))
    #er = calculateNormFullDists(er, n)
    #    cf[i] = (n > 20000) ? (n/appxInvTrace(L2;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    return distances
end

function erINV(G)
    n = G.n;
    L = lapl(G)
    Linv = mppinv(L)
    #er = diag(Linv, 0)
    min = Linv[1,1]
    u = 1
    for i = 2:n
        if Linv[i,i] < min
            u = i
            min = Linv[i,i]
        end
    end
    L2 = delnode2(L,u,n)
    Linv = inv(L2)
    cf = n/trace(Linv)
    return u, cf;
end

function erINV(G, alldistances)    
    n = G.n;
    er = zeros(n)
    L = lapl(G)
    #Linv = mppinv(L)
    u = 1
    L2 = delnode2(L,u,n)
    # inv or mppinv? TODO!
    Linv = inv(L2)
    distances = calculateCommuteDists(Linv, n, u)
    return distances
end

