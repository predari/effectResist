include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("bridge.jl")
include("components.jl")
include("compDistances.jl")
using Laplacians
using LightGraphs
#using LightGraphs.SimpleEdge

struct Component
    A:: SparseMatrixCSC{Float64}
    nc::Int64
    nodemap::Array{Int64,1} # size of nc
    bdry::Array{Int64,1} # bdry nodes in Cluster
    bdryc::Int64
end

# struct SimpleEdge{T<:Integer} #<: AbstractSimpleEdge{T}
#      src::T
#      dst::T
# end

struct Bridges
#    edges :: Array{SimpleEdge{Int64},1} #bridges
    edges:: Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}
    m :: Int64 # number of edges in Bridges
    nodes :: Set{Int64} # nodes related to Bridges instead of Array{Int64,1}
    n :: Int64 # number of nodes
    comp ::Array{Int64,1} # size of n , corresponding component
end

Component(A::SparseMatrixCSC{Float64},nodemap::Array{Int64,1}) = Component(A, A.n, nodemap,
                                                                           nothing, 0)
Bridges(edges, nodes) =
    Bridges(edges, size(edges,1),
            nodes,
            size(nodes,1),
            zeros(Int64,size(nodes,1))
            )

Bridges(edges, nodes:: Set{Int64}) =
    Bridges(edges, size(edges,1),
            nodes,length(nodes),
            zeros(Int64,length(nodes))
            )

function printComponent(C:: Component)
    println("- A type:",eltype(C.A))
    println(C.A)
    println("- A.n=",C.nc)
    println("- list of nodes=",C.nodemap)
    println("- list of bdry=",C.bdry)
    println("- bdry.n=",C.bdryc," (",100*C.bdryc/C.nc,"%)")
end

function printEdges(edges:: Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1})
    m = size(edges,1)
    
    print("- edges: ")
    for i = 1:m
        print("(",edges[i].src,",",edges[i].dst,") ")
    end
    println("")
end


function printBridges(B:: Bridges)
    #println()
    printEdges(B.edges)
    println("- m=",B.m)
    println("- list of nodes=",B.nodes)
    println("- n=",B.n)
    println("- comp of each node=",B.comp)
end

function delextnode(a::Array{Float64}, node::Int64)
    a2 = filter!(e->e!=node,a)
    return a2
end

function delextnode(a::Array{Float64}, node:: Array{Int64,1}, len :: Int64)
    for i in 1:len
        a = filter!(e->e!=node[i],a)
    end
    return a
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




function LinvDistance(a::SparseMatrixCSC{Float64}, bdry::Array{Int64,1}, bdryc::Int64 ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    if 0 in bdryc
        println("Error: node index cannot be zero!");
        exit(0)
    end
    f = approxCholLap(a,tol=1e-5);
    n = size(a,1)
    k = round(Int, JLfac*log(n)) # number of dims for JL    
    U = wtedEdgeVertexMat(a)
    m = size(U,1)

    er = zeros(n)
    er2 = zeros(n)
    cf = zeros(n, bdryc + 1)
    ## first position for d(u,v not in bdry). For d(u,bdry)
    ## position in cf is the same as in the bdry array

    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(n)
        v = f(ur[:])
        er.+= v./k
        er2.+= v.^2/k
    end

    sumer2 = sum(er2)
    sumer = sum(er)

    for j in 1:bdry
        cf[bdry[j],1] = sumer2 + n*er2[bdry[j]] -2er[bdry[j]]*sumer
        for i in 1:n
            cf[i,j + 1] = er2[i] + er2[bdry[j]] -2er[i]*er[bdry[j]]
        end
    end
    
    sumdeler2 = sum(delextnode(er2,bdry,len))
    sumdeler = sum(delextnode(er,bdry),len)
    for i in 1:n
        if (i in bdry) == false
            cf[i,1] = sumdeler2 + (n-1)*er2[i] -2er[i]*sumdeler
        end
    end
    return cf
end


function localApprox(A, brg :: Bridges, w :: IOStream)
    logw(w,"****** Running approx ******")
    n = A.n
    distances = LinvDistance(A,brg.nodes,brg.n;JLfac=200)
    return sum(distances,2);
end

function localApprox(A, extnode, w :: IOStream)
    logw(w,"****** Running approx ******")
    n = A.n
    distances = LinvDistance(A,extnode;JLfac=200)
    return sum(distances,2);
end

function removeBridges(A :: SparseMatrixCSC{Float64}, brs, nbrs)
    nodes =  Array{Int64, 1}()
    for e in brs
        A[e.src,e.dst] = 0.0
        A[e.dst,e.src] = 0.0
        push!(nodes,e.src)
        push!(nodes,e.dst)
    end
    return dropzeros!(A), Set(nodes)
end


function localApproxTest(G, bridges, w :: IOStream)
    A =  sparseAdja2(G)
    return calculateCF(localApprox(A, 0, w), A.n)
end

function cfcAccelerate(G, w :: IOStream)
    A =  sparseAdja2(G)
    n = G.n
    edges = bridges(LightGraphs.Graph(A))
    #println(brs)
    nedges = size(edges,1)
    #println(nbrs)
    #println(100*nbrs/n, "% edges")
    A, extnodes = removeBridges(A, edges, nedges)
    B = Bridges(edges, extnodes)
    cmps, map, ncmps = allComp(A)
    C = Array{Component,1}(ncmps)
    maxc = 0;
    for i in 1:ncmps
        bdry = intersect(map[i], extnodes)
        index = findin(B.nodes,bdry)
        B.comp[index] = i 
        C[i] = Component(cmps[i],cmps[i].n,map[i],bdry,size(bdry,1))
        if cmps[i].n >= maxc
            maxc = i
        end
    end
    if B.m + 1 != ncmps
        println("Error!")
        exit(0)
    end
    println("** Rate of bdrynodes compared to all nodes:", B.n*100/n,"% (",B.n,"/",n,")")
    println("** Number of components:", ncmps)
    println("** Max component size:", C[maxc].nc," bdrynodes:", C[maxc].bdryc/C[maxc].nc,"%")
    #println("Bridges:")
    #printBridges(B)
    #println("Components[",ncmps,"]:")
    #for i in 1:ncmps
    #       printComponent(C[i])
    #end
    #    
    # sizes = map(x->length(x), mapping)
    # sl = sortperm(sizes)
    # println(sizes)
    # println(sl)

    for c in C
        
    end
    
    
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
    #cf = localApproxTest(G, 0, w)
    #logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
    cfcAccelerate(G, w)
end

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
