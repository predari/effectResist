include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("bridge.jl")
include("components.jl")
include("compDistances.jl")
using Laplacians
using LightGraphs

function delextnode(a::Array{Float64}, node)
    a2 = filter!(e->e!=node,a)
    return a2
end

function LinvDistance(a::SparseMatrixCSC{Float64}, extnode::Integer; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    f = approxCholLap(a,tol=1e-5);
    n = size(a,1)
    k = round(Int, JLfac*log(n)) # number of dims for JL    
    U = wtedEdgeVertexMat(a)
    m = size(U,1)

    er = zeros(n)
    er2 = zeros(n)
    cf = zeros(n,2)

    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(n)
        v = f(ur[:])
        er.+= v./k
        er2.+= v.^2/k
    end
    if extnode != 0
        for i in 1:n
            cf[i,2] = er2[i] + er2[extnode] -2er[i]*er[extnode]
        end
    end
    sumer2 = sum(delextnode(er2,extnode))
    sumer = sum(delextnode(er,extnode))
    for i in 1:n
        cf[i,1] = sumer2 + (n-1)*er2[i] -2er[i]*sumer
    end
    return cf
end


function localApprox(A, bridges, w :: IOStream)
    logw(w,"****** Running approx ******")
    n = A.n
    distances = LinvDistance(A,0,;JLfac=200)
    return sum(distances,2);
end

function removeBridges(A :: SparseMatrixCSC{Float64}, bridges)
    for e in bridges
        A[e.src,e.dst] = 0.0
        A[e.dst,e.src] = 0.0
    end
    return dropzeros!(A)
end


function localApproxTest(G, bridges, w :: IOStream)
    A =  sparseAdja2(G)
    return calculateCF(localApprox(A, 0, w), A.n)
end

function cfcAccelerate(G, w :: IOStream)
    A =  sparseAdja2(G)
    n = G.n
    brs = bridges(LightGraphs.Graph(A))
    #println(brs)
    nbrs = size(brs,1)
    #println(nbrs)
    #println(100*nbrs/n, "% edges")
    A = removeBridges(A,brs)
    cmps, mapping, ncmps = allComp(A)
    
    if nbrs + 1 != ncmps
        println("Error!")
        exit(0)
    end
    sizes = map(x->length(x), mapping)
    sl = sortperm(sizes)
    println(sizes)
    println(sl)
    #println(cmps)
    #println(mapping)
    #println(ncmps)
    # ncomps = 1
    # println("Running components:")
    # for c in comps
    #     println(c,ncomps)
    #     if c.n != 1
    #         C = sparsemat2Graph(c)
    #         idx = nodes[ncomps]
    #         println(idx)
    #         ldistance = exact(C,1,w)
    #         ldistance = full(ldistance)
    #         s = size(ldistance,1)
    #         println("component size:",s, ",", size(ldistance,2))
    #         println(ldistance)
    #         for i in 1:s
    #             distances[idx[i],idx[i+1:end]] = ldistance[i]
    #             println(idx[i],",",idx[i+1:end])
    #             println(ldistance[i])
    #         end
    #     end
    #     ncomps = ncomps + 1
    # end
    # for i in 1:n
    #     for j in 1:n
    #         if distances[i,j] != 0.0
    #             distances[j,i] = distances[i,j]
    #         end
    #     end
    # end

    
end

datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

# for rFile in filter( x->!startswith(x, "."), readdir(string(datadir)))
#     logw(w, "---------------------",rFile,"-------------------------- ")
#     logw(w, "Reading graph from edges list file ", rFile)
#     G = read_file(string(datadir,rFile))
#     logw(w, "\t G.n: ", G.n, "\t G.m: ", G.m)
#     A, L = sparseAdja(G)
#     if !Laplacians.isConnected(A)
#         logw(w," WARNING: Graph is not connected. Program will exit!");
#         exit()
#     end
#     cf = localApproxTest(G, 0, w)
#     logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
#     cfcAccelerate(G, w)
# end

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1
logw(w, "-------------------------------------------------------- ")
n = 8
m = 20
Line = Components3Graph(n,m)
println(Line)
cf = localApproxTest(Line,0, w)
logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
println(cf)
cfcAccelerate(Line,w)
#distances = cfcaccelerate(Line, w)
#println(distances)
logw(w, "-------------------------------------------------------- ")



# function cfcAccelerate(G, w :: IOStream)
#     A =  sparseAdja2(G)
#     n = G.n
#     br = bridges(LightGraphs.Graph(A))
#     println(br)
#     #distances = zeros(n,n)
#     # distances = Array{Array{Float64, 1}}(n-1)
#     # for i in indices(distances,1) distances[i] = [] end
#     # for i in 1:n-1
#     #     for j in i+1:n
#     #         push!(distances[i], 0.0)
#     #     end
#     # end
#     A = removeBridges(A,br)
#     # for e in bg
#     #     x = e.src
#     #     y = e.dst
#     #     A[x,y] = 0
#     #     A[y,x] = 0
#     # end
#     # dropzeros!(A)

#     cmps, mapping, ncmps = allComp(A)
#     println(cmps)
#     println(mapping)
#     println(ncmps)
#     ncomps = 1
#     println("Running components:")
#     for c in comps
#         println(c,ncomps)
#         if c.n != 1
#             C = sparsemat2Graph(c)
#             idx = nodes[ncomps]
#             println(idx)
#             ldistance = exact(C,1,w)
#             ldistance = full(ldistance)
#             s = size(ldistance,1)
#             println("component size:",s, ",", size(ldistance,2))
#             println(ldistance)
#             for i in 1:s
#                 distances[idx[i],idx[i+1:end]] = ldistance[i]
#                 println(idx[i],",",idx[i+1:end])
#                 println(ldistance[i])
#             end
#         end
#         ncomps = ncomps + 1
#     end
#     for i in 1:n
#         for j in 1:n
#             if distances[i,j] != 0.0
#                 distances[j,i] = distances[i,j]
#             end
#         end
#     end
    
#     println(ncomps)
#     idx = []
#     for e in bg
#         x = e.src
#         y = e.dst
#         println(x,",",y)
#         #pull!(nbr[x][y])
#         if x < y 
#             distances[x,y] = 1
#             distances[y,x] = 1
#             for i in 1:ncomps
#                 j = find(e -> e == x ,nodes[i])
#                 if isempty(j) == false
                    
#                     c = comps[i]
#                     idx = nodes[i]
#                     #println("(",c,",",idx,")")
#                     break; # (i,j) (component, idx in nodes)
#                 end
#             end
#             #println("(",i,",",j,")")
#             idxr = filter!(e->e!=x,idx)
#             println(idxr)
#             for l in idxr
#                 distances[l,y] = distances[l,x] + 1
#             end
#             for i in 1:n
#                 for j in 1:n
#                     if distances[i,j] != 0.0
#                         distances[j,i] = distances[i,j]
#                     end
#                 end
#             end
#         end
#     end
#     distances
# return distances
# end
