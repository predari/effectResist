include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("bridge.jl")
include("components.jl")
include("compDistances.jl")
include("appxInvTrace.jl")
include("Linvdiag.jl")

include("LinvDistance.jl")
include("ShortestPaths.jl")
include("structures.jl")
using Laplacians
using LightGraphs
#using LightGraphs.SimpleEdge
using DataStructures

function isTree(A::SparseMatrixCSC)
    isConnected(A) && (nnz(A) == 2*(A.n-1))
end

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


function erJLT(A:: SparseMatrixCSC{Float64},L:: SparseMatrixCSC{Float64})

    n = A.n
    start_approx_time = time()        

    er = LinvdiagSS(A;JLfac=20)
    u = indmin(er)
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)
    cf = (n > 2000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    end_approx_time = time()
    println("TOTAL APPROX TIME ISSSS: ", end_approx_time - start_approx_time, "(s)")    
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
    #println(distances)
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
    return distances
end


function delnodes(A:: SparseMatrixCSC{Float64}, v::Set{Int64})
    idx = setdiff(Array(1:A.n),v)
    return A[idx, idx]
end
# TODO: remove
function delnode2(L, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end


function removeBridges(A :: SparseMatrixCSC{Float64}, brs, nbrs :: Integer)
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

function removeBridges(A :: SparseMatrixCSC{Float64}, B :: Bridges, core1nodes :: Set{Int64})
    ### delnodes creates a new array so the numbering
    ### is reevaluated. I want the numbering to stay the
    ### same, so I will just write zeros ontop of A whereever is needed.
    ### A = delnodes(A, B.core1nodes)
    #println(full(A))
    j = Int64
    rows = rowvals(A)
    ## or   colptr = mat.colptr , rowval = mat.rowval
    #vals = nonzeros(A)
    for u in core1nodes
        j = nzrange(A, u)[1]
        A[u,rows[j]] = 0.0
        A[rows[j],u] = 0.0
    end
    # for e in B.edges        
    #     A[e.src,e.dst] = 0.0
    #     A[e.dst,e.src] = 0.0
    # end
    for i in 1:2:length(B.edges)
        A[B.edges[i],B.edges[i+1]] = 0.0
        A[B.edges[i+1],B.edges[i]] = 0.0
    end
    return dropzeros!(A)
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

function extractBridges(A :: SparseMatrixCSC{Float64})
    start_time = time()
    B = Bridges
    B, core1nodes = bridges(LightGraphs.Graph(A))
    #println((100*B.m)/(nnz(A)/2), "% edges are bridges (type core2).")
    #println(100*length(core1nodes)/(nnz(A)/2), "% edges are bridges (type core1).")
    println((100*length(B.edges))/A.n, "% nodes are core 2++.")
    println(100*length(core1nodes)/(A.n), "% nodes are core1.")
    println("finding bridges time: ", time() - start_time, "(s)")
    #println("core1nodes",core1nodes)
    #println("Bridges:")
    #printBridges(B)
    #println("- list of core1nodes=", core1nodes)
    t = time()
    A  = removeBridges(A, B, core1nodes)
    println("remove bridges time: ", time()- t, "(s)")
    return A,B
end

function buildComponents(A :: SparseMatrixCSC{Float64}, B :: Bridges)
    start_time = time()
    cmps, map, ncmps = allComp(A)
    println("finding components time: ", time() - start_time, "(s)")
    t = time()
    C = Array{Component,1}(ncmps)
    maxc = 0;
    #### Is : for i = eachindex(a) faster than for i = 1:n?

    for i in eachindex(cmps) #1:ncmps
        bdry = collect(intersect(DataStructures.SortedSet(B.core2nodes), Set(map[i])))
        # TODO: check if link has to be ordered or not!!!!
        # TODO: check if intersect is still that bad with (sorted) arrays (compared to Sets)
        #link = collect(intersect(DataStructures.SortedSet(B.edges), Set(map[i])))
        link = collect(intersect(Set(B.edges), Set(map[i])))
        #println("in buildComp")
        #println(bdry)
        #println(link)

        index = findin(B.core2nodes,bdry)
        #println(B.ext[index])
        B.comp[findin(B.edges,link)] = i 
        C[i] = Component(cmps[i],cmps[i].n,map[i],bdry,link,length(link),
                         zeros(cmps[i].n,length(link)), B.ext[index])
        if cmps[i].n >= maxc
            maxc = i
        end
    end
    println("building structure components time: ", time() - t, "(s)")
    return C
end

function stripcore1nodes(comp  :: Array{Int64,1} , core3nodes :: Array{Int64,1})
    ## assert such that length(comp) and length(core3nodes) is equal
    pst = []
    for i in eachindex(comp)
        if comp[i] == 0
            push!(pst,i)
        end
    end
    #println(pst)
    deleteat!(comp,pst)
    deleteat!(core3nodes,pst)
    return comp, core3nodes
end
function createTreeGraph(n  :: Int64)
        # the construction of this graph is buggy!
        ### TODO: think about what kind of graph (chinese type, sparsematrixCSC, lightgraph etc)
        ### graph type you need!!
        nbr = Array{Array{Int, 1}}(n)
        for i in indices(nbr,1) nbr[i] = [] end
        wgt = Array{Array{Int, 1}}(n)
        for i in indices(wgt,1) wgt[i] = [] end
        for i in 1:n-1
            push!(nbr[i], i + 1)
            push!(nbr[i + 1],i)
            # if ncomp[i] != ncomp[i + 1]
            #     wgt[] = something
            # else
            #     wgt[] = 1
            # end
        end
    return Graph(n,n-1,nbr)
end



function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream)
    start_time = time()
    n = A.n
    B = Bridges 
    A, B = extractBridges(A)
    # for (i,e) in enumerate(B.edges)
    #     B.core3nodes[2*i - 1] = e.src
    #     B.core3nodes[2*i] = e.dst
    # end
    t = time()
    C = Array{Component,1}
    C = buildComponents(A, B)
    count = length(C)
    println("creating components time: ", time()- t, "(s)")

    t = time()
    
    for (idx, c) in enumerate(C)
        print("Approxing component $idx ...")
        c.distances = localApprox(c, w)
        println(" done")
        #cf = exact(sparsemat2Graph(c.A), w )
    end
    println("Bridges:")
    printBridges(B)
    println("Components: $count")
    for (idx, c) in enumerate(C)
        print("$idx")
        printComponent(c)
    end

    println(" solving core1 time : ", time()- t, "(s)")
    if count == 1
        c = C[1]
        c.distances = sum(c.distances,2)
        cf = calculateCF(c.distances, n, c.nc)
        #println("Final", cf)
        logw(w,"\t node with argmax{c(", c.nodemap[indmax(cf)], ")} = ", maximum(cf))
    else
        # this is not correct yet, need to also remove the info for nodes
        # connected with core1nodes
        println("length = ",length(B.edges))
        println("components:", B.comp)
        println("edges:", B.edges)
        
        newcomp, newedges = stripcore1nodes(B.comp,B.edges)
        println(newcomp)
        println(newedges)
        n = length(newcomp)
        # cG = createTreeGraph(n)
        # println(cG)
        # cA, cL = sparseAdja(cG)
        cVertices = Array{cVertex}(count)
        nodes = Array{Array{Int, 1}}(count)
        rnodes = Array{Array{Int, 1}}(count)
        for i in indices(nodes,1) nodes[i] = [] end
        for i in indices(rnodes,1) rnodes[i] = [] end
        for (i,e) in enumerate(newedges)
            push!(rnodes[newcomp[i]], e)
            push!(nodes[newcomp[i]], i)
        end
        #println("cVertices:")
        for i in 1:count
            cVertices[i] = cVertex(i,nodes[i],rnodes[i],count)
            #printcVertex(cVertices[i])
        end
        cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
        for i in 1:n-1
            if newcomp[i] != newcomp[i + 1]
                cA[i,i+1] .= 1.0
                cA[i+1,i] .= 1.0
            else
                c = C[newcomp[i]]

                #println(c.link, newedges[i])
                lidx1 = findin(c.link, newedges[i])
                lidx2 = findin(c.nodemap, newedges[i + 1])
                #### TODO: restore values!!!
                dist = getindex(c.distances[lidx1 + 1,lidx2])
                #println("idx = ",lidx1 + 1,lidx2, dist)
                cA[i,i+1] .= dist
                cA[i+1,i] .= dist
            end
        end
        if !isTree(cA)
            logw(w,"WARNING: contracted graph should be a tree! ");
        end
        #println(cA)
        println("number of components:", n)
        
        dist2 = zeros(Float64,count,count) # or count
        for i in 1:n
            x = newcomp[i]
            if sum(dist2[x,:]) != 0.0
                tmp , p =  shortestPaths(cA, i, newcomp, cVertices, count)
                #shortestPaths(cA, i, newcomp, count)
                println(p)
                for j in 1:count
                    dist2[x,j] = dist2[x,j] < tmp[j] ? dist2[x,j] : tmp[j] 
                end
            else
                dist2[x,:], p = shortestPaths(cA, i, newcomp, cVertices, count)
                #shortestPaths(cA, i, newcomp, count)
            end
            #println("dist2",dist2)
            #println("p",p)
        end
        println("Final dist2:")
        println(dist2)
        println("dist2 after summation:")
        distcomp = sum(dist2,2)
        println(distcomp)
        rx :: Int64 = 0
        for i in 1:count
            # for i I need Component structure (C[])
             for j in 1:count
                # for j I need cVertex structure (cVertices[])
                if j == i continue; end
                if cVertices[i].cmplist[j] == 0
                    rx = getindex(cVertices[j].rx[1])
                else
                    idx = findin(cVertices[j].x, cVertices[i].cmplist[j])
                    rx = getindex(cVertices[j].rx[idx])
                    #rx = cVertices[i].cmplist[j]
                end
                idx = findin(C[j].link, rx)
                #println("rx and idx", rx, idx+1)
                #println(C[j].distances[:,idx+1], sum(C[j].distances[:,idx+1]))
                distcomp[i] += sum(C[j].distances[:,idx+1])
             end
        end
        println("distcomp:")
        println(distcomp)
        for (idx, c) in enumerate(C)
            if length(c.link) == 2
                for i in length(c.link)
                    ## add it one more time
                    c.distances[:,i+1] += c.distances[:,i+1]
                end
            elseif length(c.link) == 1
                ## add it two more times (this is hardcoded, I have to find how many times)
                c.distances[:,2] += 2*c.distances[:,2]
            end
            c.distances = sum(c.distances,2) + distcomp[idx]
            println("distances:")
            println(c.distances)
            #c.distances = sum(c.distances,2)
            #cf = calculateCF(c.distances, n, c.nc)
        end
    end
    
    # # ncomps = 1
    # # for c in comps
    # #     if c.n != 1
    # #         C = sparsemat2Graph(c)
    # #         idx = nodes[ncomps]
    # #         ldistance = exact(C,1,w)
    # #         ldistance = full(ldistance)
    # #         s = size(ldistance,1)
    # #         for i in 1:s
    # #             distances[idx[i],idx[i+1:end]] = ldistance[i]
    # #             println(idx[i],",",idx[i+1:end])
    # #         end
    # #     end
    # #     ncomps = ncomps + 1
    # # end
    println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
end

datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1
logw(w, "-------------------------------------------------------- ")
#n = 8
#m = 20
#Line = Components3Graph(n,m)

Line = TestGraph(18, 48)
println(Line)
A, L = sparseAdja(Line)
@time approx(A,L,w)
Line = TestGraph(18, 48)
println(Line)
A, L = sparseAdja(Line)
@time cfcAccelerate(A,w)
@time exact(Line,w)

# Line = subTestGraph(11, 30)
# println(Line)
# A, L = sparseAdja(Line)
# cfcAccelerate(A,w)
# approx(A,L,w)
# Line = subTestGraph2(8, 24)
# println(Line)
# A, L = sparseAdja(Line)
# cfcAccelerate(A,w)
# #println("TODO: here you multiply by: ",Line.n,"instead of the size of cf: ", size(cf,1))
# #cf = calculateCF(cf, Line.n, size(cf,1))
# #logw(w,"\t node with argmax{c(", indmax(cf), ")} = ", maximum(cf))
# #distances = cfcaccelerate(Line, w)
# #println(distances)
# approx(A,L,w)
# logw(w, "-------------------------------------------------------- ")

exit()

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
    @time cfcAccelerate(A, w)
    A, L = sparseAdja(G)
    @time approx(A,L,w)
    #@time exact(G,w)
    exit()
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
