include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("bridge.jl")
include("components.jl")
include("compDistances.jl")
include("appxInvTrace.jl")
include("Linvdiag.jl")

include("LinvDistance.jl")
include("structures.jl")
using Laplacians
using LightGraphs
#using LightGraphs.SimpleEdge

function delextnode(a::Array{Float64}, node::Int64)
    a2 = filter!(e->e!=node,a)
    return a2
end

function delextnode(a::Array{Float64}, node:: Array{Int64,1}, len :: Int64)
    a2 = a
    for i in 1:len
        a2 = filter!(e->e!=node[i],a2)
    end
    return a2
end

function checkDistancesComponent(cf, A :: SparseMatrixCSC{Float64})
    println(cf)
    res1 = calculateCF(cf, A.n)
    cf2 = localApprox(A, 0, w)
    println(cf2)
    res2 = calculateCF(cf2, A.n)
    if res1 == res2
        println(":) - The distance calculations are correct!")
    else
        println(":( - Please check again your distances calculations!") 
    end
end



# TODO: remove
function delnode2(L, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function erJLT(A:: SparseMatrixCSC{Float64},L:: SparseMatrixCSC{Float64})

    n = A.n
    start_time = time()        
    er = LinvdiagSS(A;JLfac=20)
    u = indmin(er)
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)
    cf = (n > 2000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    
    println("TOTAL APPROX TIME IS: ", time() - start_time, "(s)")    
  return u, cf;

end

function approx(A:: SparseMatrixCSC{Float64},L:: SparseMatrixCSC{Float64}, w :: IOStream)
    logw(w,"****** Running (chinese) approx ******")
    u, maxcf = erJLT(A,L)
    logw(w,"\t node with argmax{c(", u, ")} = ", maxcf)
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
    println(distances)
    return distances
end


function exact(G, w :: IOStream)
    logw(w,"****** Running (my) exact ******")
    distances = erINV(G,1)
    cf = calculateNodeDists(distances, G.n)
    println("Exact distance:",cf)
    cf = calculateCF(cf, G.n)
    logw(w,"\t node with argmax{c(", indmax(cf), ")} = ",
         maximum(cf))
    return cf
end

#function localApprox(A, brg :: Bridges, w :: IOStream)
function localApprox(A, extnodes:: Array{Int64,1}, size::Integer , w :: IOStream)
    logw(w,"****** Running approx ******")
    n = A.n
    distances = LinvDistance(A,extnodes,size;JLfac=200)
    #println(distances)
    return sum(distances,2);
end

function localApprox(A, extnode:: Integer, w :: IOStream)
    logw(w,"****** Running approx ******")
    n = A.n
    distances = LinvDistance(A,extnode;JLfac=200)
    #println(distances)
    return sum(distances,2);
end

function localApprox(c :: Component, w :: IOStream)
    distances = LinvDistance(c)
    #distances = LinvDistance(c.A, c.bdry, c.bdryc, c.nodemap)
    #println(distances)
    distances = sum(distances,2)
    #println(distances)
    c.distances = distances
    return distances
end

function removeBridges(A :: SparseMatrixCSC{Float64}, brs, nbrs)
    #nodes =  Array{Int64,1}()
    nodes =  Set{Int64}()
    for e in brs
        A[e.src,e.dst] = 0.0
        A[e.dst,e.src] = 0.0
        push!(nodes,e.src)
        push!(nodes,e.dst)
    end
    return dropzeros!(A), nodes
end


function localApproxTest(G, bridges, w :: IOStream)
    A =  sparseAdja2(G)
    return calculateCF(localApprox(A, 0, w), A.n)
end

function locateBridges(A :: SparseMatrixCSC{Float64})
    edges = bridges(LightGraphs.Graph(A))
    #println("Time of finding bridges!")
    #nedges = size(edges,1)
    #println((100*nedges)/(G.m), "% edges are bridges.")
    A, extnodes = removeBridges(A, edges, nedges)
    B = Bridges(edges, extnodes)
    return A,B
end

function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream)
    start_time = time()
    n = A.n
    edges = bridges(LightGraphs.Graph(A))
    nedges = size(edges,1)
    println((100*nedges)/(A.m), "% edges are bridges.")
    println("finding bridges time: ", time() - start_time, "(s)")
    t = time()
    A, extnodes = removeBridges(A, edges, nedges)
    B = Bridges(edges, extnodes)
    cmps, map, ncmps = allComp(A)
    C = Array{Component,1}(ncmps)
    println("remove bridges time: ", time()- t, "(s)")
    t = time()
    maxc = 0;
    for i in 1:ncmps
        # Intersect with arrays is slow because in is slow with arrays.
        # intersect(Set(a),Set(b)) is better. Or setdiff(a, setdiff(a, b))
        bdry = collect(intersect(Set(map[i]), extnodes)) 
        index = findin(B.nodes,bdry)
        B.comp[index] = i 
        C[i] = Component(cmps[i],cmps[i].n,map[i],bdry,length(bdry),
                         zeros(cmps[i].n,length(bdry)), zeros(length(bdry)))
        if cmps[i].n >= maxc
            maxc = i
        end
    end
    println("creating components time: ", time()- t, "(s)")
    # if B.m + 1 != ncmps
    #     println("Error!")
    #     exit(0)
    # end
    #### stage one now
    #### for all bridges that connect one-core components with other components
    t = time()
    for (idxu, u) in enumerate(B.nodes)
        v = u
        i = B.comp[idxu]
        if C[i].nc == 1
            for e in B.edges
                if e.src == u
                    v = e.dst
                    break;
                elseif e.dst == u
                    v = e.src
                    break;
                end
                #todo:remove e with pop!
            end
            tmp = C[B.comp[findin(B.nodes,v)]]
            c=tmp[1]
            c.external[findin(c.bdry, v)] = c.external[findin(c.bdry, v)] + 1
        end 
    end
    println(" creating core1 time : ", time()- t, "(s)")
    ### end of stage one
    t = time()
    for (idx, c) in enumerate(C)
        if c.nc != 1
            #cf = localApprox(c.A, 0, w)
            print("Approxing component $idx ...")
            cf = localApprox(c, w)
            println(" done")
            #cf = exact(sparsemat2Graph(c.A), w )
        end
    end
    println(" solving core1 time : ", time()- t, "(s)")
    # for (idx, c) in enumerate(C)
    #     if c.nc == 1
    #         println("$idx")
    #         deleteat!(C,idx)
    #     end
    # end
    
    # ncomps = 1
    # for c in comps
    #     if c.n != 1
    #         C = sparsemat2Graph(c)
    #         idx = nodes[ncomps]
    #         ldistance = exact(C,1,w)
    #         ldistance = full(ldistance)
    #         s = size(ldistance,1)
    #         for i in 1:s
    #             distances[idx[i],idx[i+1:end]] = ldistance[i]
    #             println(idx[i],",",idx[i+1:end])
    #         end
    #     end
    #     ncomps = ncomps + 1
    # end
    println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
end

datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1
logw(w, "-------------------------------------------------------- ")
# # n = 8
# # m = 20
# # Line = Components3Graph(n,m)
# # Line = TestGraph(15, 40)
# # println(Line)
# # @time cf = localApproxTest(Line,0, w)
# # #println(cf)
# # logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
# # @time cfcAccelerate(Line,w)
# Line = subTestGraph(11, 30)
# println(Line)
# @time cfcAccelerate(Line,w)
# #println("TODO: here you multiply by: ",Line.n,"instead of the size of cf: ", size(cf,1))
# #cf = calculateCF(cf, Line.n,size(cf,1))
# #logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
# @time cf = localApproxTest(Line,0, w)
# logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))

# Line = subTestGraph2(8, 24)
# println(Line)
# @time cfcAccelerate(Line,w)
# #println("TODO: here you multiply by: ",Line.n,"instead of the size of cf: ", size(cf,1))
# #cf = calculateCF(cf, Line.n, size(cf,1))
# #logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
# #distances = cfcaccelerate(Line, w)
# cf = localApproxTest(Line,0, w)
# logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
# #println(distances)
# logw(w, "-------------------------------------------------------- ")


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

    cfcAccelerate(A, w)
   
    approx(A,L,w)

    #@time cf = localApproxTest(G, 0, w)
    #logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
end




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
