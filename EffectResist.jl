include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("bridge.jl")
include("components.jl")
include("Exact.jl")
include("appxInvTrace.jl")
include("Linvdiag.jl")
using Laplacians
using LightGraphs

function isTree(A::SparseMatrixCSC)
    isConnected(A) && (nnz(A) == 2*(A.n-1))
end


function erJLT(G)
    n = G.n
    
    A, L = sparseAdja(G)
    er = LinvdiagSS(A;JLfac=20)
    u = indmin(er)
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)
    cf = (n > 20000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
  return u, cf;

end


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
        println(typeof(V))
        println(size(V,1),",",size(V,2))        
        #println(typeof(r))
        #V[:,i] .= V[:,i]*R[:,i]
        #rst += sum(v.*r)
    end
    return V/k
end

function erJLT(G, alldistances)
    n = G.n
    cf = zeros(n)
    A, L = sparseAdja(G)
    #er = LinvdiagSS(A;JLfac=20)
    #u = indmin(er)
    u = 1
    #L2 = delnode2(L,u,n)
    A2 = delnode2(A,u,n)
    Linv = LinvSS(A2;JLfac=20)
    # TODO: distances should be stored as adj list.
    println(typeof(Linv))
    distances = calculateEdgeDists(Linv, n, u)
    #    cf[i] = (n > 20000) ? (n/appxInvTrace(L2;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    return distances
end

function approx(G, alldistances, w :: IOStream)
    logw(w,"****** Running approx ******")
    distances = erJLT(G, 1)
    if alldistances == 1
        return distances
    else
        cf = calculateNodeDists(distances, G.n)
        println(cf)
        cf = calculateCF(cf, G.n)
        logw(w,"\t node with argmax{c(", indmax(cf), ")} = ",
             maximum(cf))
        return cf
    end
end

function approx(G, w :: IOStream)
    logw(w,"****** Running (chinese) approx ******")
    #start_time = time()
    u, maxcf = erJLT(G)
    #elapsed_time = time()-start_time
    logw(w,"\t node with argmax{c(", u, ")} = ", maxcf)
    #logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")
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
    distances = calculateEdgeDists(Linv, n, u)
    return distances
end



function calculateEdgeDists(Linv, n, u) # n is the size of G not L. size(L,1) = n-1
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
    end
    return cf
end

function exact(G, alldistances, w :: IOStream)
    logw(w,"****** Running (my) exact ******")
    # time start
    # start_time = time()
    distances = erINV(G,1)
    if alldistances == 1
        # elapsed_time = time()-start_time
        # logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")
        return distances
    else
        cf = calculateNodeDists(distances, G.n)
        cf = calculateCF(cf, G.n)
        logw(w,"\t node with argmax{c(", indmax(cf), ")} = ",
             maximum(cf))
        #elapsed_time = time()-start_time
        #logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")
        return cf
    end
end

function exact(G, w :: IOStream)
    logw(w,"****** Running (chinese) exact ******")
    # time start
    #start_time = time()
    u, maxcf = erINV(G)
    #elapsed_time = time()-start_time
    logw(w,"\t node with argmax{c(", u, ")} = ", maxcf)
    #logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")

end



datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1



for rFile in filter( x->!startswith(x, "."), readdir(string(datadir)))
    logw(w, "---------------------",rFile,"-------------------------- ")
    logw(w, "Reading graph from edges list file ", rFile)
    G = read_file(string(datadir,rFile))
    logw(w, "\t G.n: ", G.n, "\t G.m: ", G.m)
    A, L = sparseAdja(G)
    if !Laplacians.isConnected(A)
        logw(w," WARNING: Graph is not connected. Program will exit!");
        exit()
    end
    @time exact(G,w)
    @time exact(G,0,w)
    @time approx(G,w)
    @time approx(G,0,w)
end
    logw(w, "-------------------------------------------------------- ")

n = 8
m = 20
Line = Components3Graph(n,m)
println(Line)
#exact_distance = exact(Line,1,w) # same as er_mppinv(Line)
exact(Line,w)
exact(Line,0,w)
    
# LightGraphs structure is needed only to locate bridges
SLine =  sparseAdja2(Line)
bg = bridges(LightGraphs.Graph(SLine))
println(bg)
dist = zeros(Float64,n,n)
for e in bg
    x = e.src
    y = e.dst
    println(x,",",y)
    #pull!(nbr[x][y])
    SLine[x,y] = 0
    SLine[y,x] = 0
    dist[x,y] = 1
    dist[y,x] = 1
end
dropzeros!(SLine)
comps, nodes = allComp(SLine)
#println(comps)
println(nodes)
println("Running components:")
for c in comps
    if c.n != 1
        # ------- exact -------
        # println(c)
        C = sparsemat2Graph(c)
        ldist = exact(C,0,w)
        # ---------------------
        # ------- approx -------
        ### C = convert(SparseMatrixCSC{Float64,Int64},c)
        # println(sparsemat2Graph(c))
        # ldist = approx(c, w)
        # ---------------------
        #println("Cf values for components")
        #println(ldist)
        
    end
end
    
# dropzero! is needed to disconnect the components for real
#     println("after deletion of bridges")
#     println(full(spPath))
#     println(spPath)
#     comps, nodes = allComp(spPath)
#     println(nodes)
#     println("components")
#     for c in comps
#         if c.n != 1
#             local_er = approx(convert(SparseMatrixCSC{Float64,Int64},c),w)
#             #global_er   
#         end
#     end

    

    #T1 = akpwU(A);
    #if !isTree(T1)
    #    logw(w,"WARNING: Graph is a tree! ");
    #end
    #comp = randperm(MersenneTwister(1234), A.n)
    #compGraphU(A, comp)
    logw(w,"****** Running starDecomp ******")
    #T = starDecomp(A)
    #print(T)
    # if !isTree(T)
    #     logw(w,"WARNING: Graph is a tree! ");
    # end
    # logw(w, "\t T.n: ", T.n, "\t T.m: ", nnz(T))

close(w)
