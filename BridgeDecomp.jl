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
    #println("Exact distance:",distances)
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
    #println(core1nodes)
    println((100*B.m)/(nnz(A)/2), "% edges are bridges type core2.")
    println(100*length(B.edges)/A.n, "% nodes are core2.")
    println(100*length(core1nodes)/(nnz(A)/2), "% edges are bridges type core1.")
    println( 100*length(core1nodes)/(A.n), "% nodes are core1.")
    println("finding bridges time: ", time() - start_time, "(s)")
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


function contractAdjGraph(edges :: Array{Int64,1}, cmplist :: Array{Int64,1} , C :: Array{Component,1}, nc :: Int64)
    n = length(cmplist)
    cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
    for i in 1:2:length(edges)
        cA[i,i+1] .= 1.0
        cA[i+1,i] .= 1.0
    end
    #println(cA)
    #nc = length(C)
    for i in 1: nc
        idx = findin(cmplist,i)
        if length(idx) > 1
            c = C[i]
            #println("component:",i)
            #println(idx, edges[idx])
            lidx1 = findin(c.link, edges[idx])
            lidx2 = findin(c.nodemap, edges[idx])
            lidx1 += 1
            #println(lidx1, lidx2)
            dist = c.distances[lidx2,lidx1]
            #println(dist)
            cA[idx,idx] = dist
        end
        #println(cA)
    end
    for i in 1:n
        if cA[i,i] != 0.0
            cA[i,i] = 0.0
        end
    end
    dropzeros!(cA)
    if !isTree(cA)
        logw(w,"WARNING: contracted graph should be a tree! ");
    end
    return cA
end

function shortestContractPaths(cA :: SparseMatrixCSC{Float64}, nc :: Int64, edges :: Array{Int64,1}, cmplist :: Array{Int64,1})
    n = length(cmplist)
    dist2 = zeros(Float64,nc,nc)
    path2 = zeros(Int64,nc,nc)
    for i in 1:n
        x = cmplist[i]
        if sum(dist2[x,:]) != 0.0
            tmp , path =  shortestPaths(cA, i, cmplist, edges, nc)
            for j in 1:nc
                dist2[x,j] = dist2[x,j] < tmp[j] ? dist2[x,j] : tmp[j] 
            end
        else
            dist2[x,:], path = shortestPaths(cA, i, cmplist, edges, nc)
        end
        if path2[x,1] == 0
            path2[x,:] = path
        end
    end
    return dist2, path2
end

function compRealSizes(C :: Array{Component,1}, nc :: Int64)
    sizes = zeros(Int64,nc)
    for i in 1:nc
        externalNodes = 0
        for j in C[i].external
            externalNodes += j 
        end
        sizes[i] = C[i].nc + externalNodes
    end
    return sizes
end
function compContractLocalDists(C :: Array{Component,1}, nc :: Int64, path2 :: Array{Int64,2}, dist2 :: Array{Float64,2}, sizes :: Array{Int64,1})
    distcomp = zeros(Float64,nc)
   
    for i in 1:nc
        println("- C[$i] cf-distances=",C[i].distances)
        # TODO: for (j,p) in enumerate(C[i].path)
        for (j,p) in enumerate(path2[i,:])
            if j == i continue; end
            idx = findin(C[j].link, p)
            println("In comp=$j node-->$idx sum:", sum(C[j].distances[:,idx+1]))
            distcomp[i] += sizes[j]*dist2[i,j] + sum(C[j].distances[:,idx+1])
        end
        #### following code for debug purposes
        # println(dist2)
        # for (j,p) in enumerate(path2[i,:])
        #     if j == i continue; end
        #     println(sizes[j],"*",dist2[i,j])
        #     distcomp[i] += sizes[j]*dist2[i,j]
        # end
        # println("distcomp after size:", distcomp)
        # for (j,p) in enumerate(path2[i,:])
        #     if j == i continue; end
        #     idx = findin(C[j].link, p)
        #     println("In comp=$j node-->$p sum:", sum(C[j].distances[:,idx+1]))
        #     distcomp[i] += sum(C[j].distances[:,idx+1])
        # end
        # println("distcomp after sum:", distcomp)
    end
    return distcomp
end


function updateLocalDists(C :: Array{Component,1}, sizes :: Array{Int64,1},  edges :: Array{Int64,1}, cmplist :: Array{Int64,1})  
    s :: Int64 = sum(sizes)
    l :: Int64 = length(cmplist)
    for (idx, c) in enumerate(C)
        if size(c.distances,2) < 2
            println("Warning: # of columns for distances should be at least 2")
        end
        if size(c.distances,2) == 2
            println("In comp=$idx link with idx=2, (mult with ", s - sizes[idx], ")")
            c.distances[:,2] += (c.distances[:,2] * (s - sizes[idx]))
        else
            for i in 2:size(c.distances,2)
                node = c.link[i-1]
                x = getindex(findin(edges,node))
                if rem.((l-x),2) == 1
                    y = x + 1
                else # (==0)
                    y = x - 1
                end
                yi:: Int64 = cmplist[y]
                println("In comp=$idx link with idx=$i, (mult with ", sizes[yi], ")")
                c.distances[:,i] += (c.distances[:,i] * sizes[yi])
            end
        end
    end
end


function aggregateLocalDists(C :: Array{Component,1}, distcomp :: Array{Float64,1})
    fdistance = []
    fnodes = []
    for (idx, c) in enumerate(C)
        c.distances = sum(c.distances,2) + distcomp[idx]
        println("distances:")
        println(c.distances)
        fdistance = [fdistance ; c.distances]
        fnodes = [fnodes ; c.nodemap]
    end
    return fdistance, fnodes
end

function aggregateDistances(C :: Array{Component,1}, nc :: Int64, path2 :: Array{Int64,2}, dist2 :: Array{Float64,2}, sizes :: Array{Int64,1}, edges :: Array{Int64,1}, cmplist :: Array{Int64,1})
    distcomp = zeros(Float64,nc)
    fdistance = []
    fnodes = []
    for i in 1:nc
        #println("- C[$i] cf-distances=",C[i].distances)
        # TODO: for (j,p) in enumerate(C[i].path)
        for (j,p) in enumerate(path2[i,:])
            if j == i continue; end
            idx = findin(C[j].link, p)
            #println(j," ",p," ",idx)
            distcomp[i] += sizes[j]*dist2[i,j] + sum(C[j].distances[:,idx+1])
        end
    end
    println("distcomp:")
    println(distcomp)
    #println(B.edges," ",length(B.edges)," ",B.comp)
    #sizes -= 1
    for i in 1:2:length(edges)
        #print("(",B.edges[i],",",B.edges[i+1],") \n")
        for (idx, c) in enumerate(cmplist[i:i+1])
            #println(idx," ", c," ", C[c].link)
            for i in 1:length(C[c].link)                   
                #println(B.comp[Int(2/idx)], ": ", sizes[B.comp[Int(2/idx)]])
                C[c].distances[:,i+1] += (C[c].distances[:,i+1] * sizes[ cmplist[Int(2/idx)]])
            end
        end
    end
    for (idx, c) in enumerate(C)
        c.distances = sum(c.distances,2) + distcomp[idx]
        println("distances:")
        println(c.distances)
        fdistance = [fdistance ; c.distances]
        fnodes = [fnodes ; c.nodemap]
    end
    return fdistance, fnodes
end




function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream)
    start_time = time()
    n = A.n
    B = Bridges 
    A, B = extractBridges(A)
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
    end
    println("Bridges:")
    printBridges(B)
    # println("Components: $count")
    # for (idx, c) in enumerate(C)
    #     print("$idx")
    #     printComponent(c)
    # end
    println(" solving core1 time : ", time()- t, "(s)")
    if count == 1
        c = C[1]
        c.distances = sum(c.distances,2)
        cf = calculateCF(c.distances, n, c.nc)
        logw(w,"\t node with argmax{c(", c.nodemap[indmax(cf)], ")} = ", maximum(cf))
    else
        println("components:", B.comp)
        println("edges:", B.edges)
        ## strip core1nodes in order to remove internal paths
        newcomp, newedges = stripcore1nodes(B.comp,B.edges)
        println("After striping nodes1!")
        println("components:", newcomp)
        println("edges:", newedges)
        n :: Int64 = length(newcomp)
        # following line: improves perf? TODO: check
        cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
        cA = contractAdjGraph(newedges, newcomp, C, count)
        println(cA)
        # following 2 lines: improves perf? TODO: check
        dist2 = zeros(Float64,count,count)
        path2 = zeros(Int64,count,count)
        dist2, path2 = shortestContractPaths(cA, count, newedges, newcomp)
        println("path2:", path2)        
        println("dist2:", dist2)
        sizes = zeros(Int64,count)
        sizes = compRealSizes(C, count)
        println("Sizes:", sizes)
        distcomp = zeros(Float64,count)
        distcomp = compContractLocalDists(C, count, path2, dist2, sizes)
        println("Final distcomp:", distcomp)
        updateLocalDists(C, sizes,newedges, newcomp)
        fdistance, fnodes = aggregateLocalDists(C, distcomp)
        #fdistance, fnodes = aggregateDistances(C, count, path2, dist2, sizes, newedges, newcomp)
        #println(fdistance)
        #println(fnodes)
        cf = calculateCF(fdistance, A.n,length(fdistance))
        logw(w,"\t node with argmax{c(", fnodes[indmax(cf)], ")} = ", maximum(cf))
        println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
    end
end

datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1
logw(w, "-------------------------------------------------------- ")
#n = 8
#m = 20
#Line = Components3Graph(n,m)

#Line = TestGraph(20, 54)
Line = TestGraph(18, 48)
#Line = Components3Graph(8, 22)
#Line = subTestGraph(14, 38)
println(Line)
A, L = sparseAdja(Line)
@time approx(A,L,w)
#Line = Components3Graph(8, 22)
Line = TestGraph(18, 48)
#Line = TestGraph(20, 54)
#Line = subTestGraph(14, 38)
println(Line)
A, L = sparseAdja(Line)
@time cfcAccelerate(A,w)
@time exact(Line,w)


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
#    @time exact(G,w)
end

