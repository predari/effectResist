
include("logw.jl")
include("LinvDistance.jl")
include("ShortestPaths.jl")
include("structures.jl")
include("bridge.jl")
include("components.jl")
include("compDistances.jl")

using Laplacians


#function localApprox(A, brg :: Bridges, w :: IOStream)
function localApprox(A, extnodes:: Array{Int64,1}, size::Integer , w :: IOStream)
    logw(w,"****** Running fast approx ******")
    n = A.n
    distances = LinvDistance(A,extnodes,size;JLfac=200)
    return sum(distances,2);
end

function localApprox(A, extnode:: Integer, w :: IOStream)
    logw(w,"****** Running fast approx ******")
    n = A.n
    distances = LinvDistance(A,extnode;JLfac=200)
    return sum(distances,2);
end

function localApprox(c :: Component, w :: IOStream)
    distances = LinvDistance(c)
    return distances
end


# function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream)
#     start_time = time()
#     n = A.n
#     #computeCore2Bridges(A)
#     B = Bridges
#     A, B = extractBridges(A)
#     t = time()
#     C = Array{Component,1}
#     C = buildComponents(A, B)
#     count = length(C)
#     for c in C
#         print(c.size,"-",length(c.nodemap)," ")
#     end
#     println("creating components time: ", time()- t, "(s)") 
#     t = time()
#     for (idx, c) in enumerate(C)
#         print("Approxing component $idx ...")
#         c.distances = localApprox(c, w)
#         println(" done")
#     end
#     # println("Bridges:")
#     # printBridges(B)
#     # println("Components: $count")
#     # for (idx, c) in enumerate(C)
#     #     print("$idx")
#     #     printComponent(c)
#     # end
#     println(" solving core1 time : ", time()- t, "(s)")
#     if count == 1
#         c = C[1]
#         c.distances = sum(c.distances,2)
#         cf = calculateCF(c.distances, n, c.nc)
#         logw(w,"\t node with argmax{c(", c.nodemap[indmax(cf)], ")} = ", maximum(cf))
#     else
#         #println("components:", B.comp)
#         #println("edges:", B.edges)
#         #println("size:",length(B.comp))
#         ## strip core1nodes in order to remove internal paths
#         ##### stripcore1nodes is in components.jl
#         newcomp, newedges = stripcore1nodes(B.comp,B.edges)
#         #println("After striping nodes1!")
#         #println("components:", newcomp)
#         #println("size:",length(newcomp))
        
#         #println("edges:", newedges)
#         n :: Int64 = length(newcomp)
#         # following line: improves perf? TODO: check
#         cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
#         cA = contractAdjGraph(newedges, newcomp, C, count)
#         #println(cA)
#         # following 2 lines: improves perf? TODO: check
#         dist2 = zeros(Float64,count,count)
#         path2 = zeros(Int64,count,count)
#         dist2, path2 = shortestContractPaths(cA, count, newedges, newcomp)
#         #println("path2:", path2)        
#         #println("dist2:", dist2)
#         sizes = zeros(Int64,count)
#         sizes = compRealSizes(C, count)
#         #println("Sizes:", sizes)
#         distcomp = zeros(Float64,count)
#         distcomp = compContractLocalDists(C, count, path2, dist2, sizes)
#         println("Final distcomp:", distcomp)
#         updateLocalDists(C, sizes,newedges, newcomp)
#         fdistance, fnodes = aggregateLocalDists(C, distcomp)
#         #fdistance, fnodes = aggregateDistances(C, count, path2, dist2, sizes, newedges, newcomp)
#         #println(fdistance)
#         #println(fnodes)
#         cf = calculateCF(fdistance, A.n,length(fdistance))
#         logw(w,"\t node with argmax{c(", fnodes[indmax(cf)], ")} = ", maximum(cf))
#         println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
#     end
# end

            

function isSpecialNode(bridges :: Array{Int64,1}, bdrynodes :: Array{Int64,1}, u = nothing)
    if u != nothing
        if in(u,bridges) == true
            println("MAX VALUE $u IS PART OF LINK NODES!!!")
        end
        if in(u,bdrynodes) == true
            println("MAX VALUE $u IS PART OF CORE2NODES NODES!!!")
        end
    end
end

### TODO: move in components.jl 
function maxComponentIdx(C :: Array{Component,1})
    max :: Int64 = 0
    idx :: Int64 = 0
    for (i, c) in enumerate(C)
        if c.size > max ### size or c.nc
            max = c.size
            idx = i
        end
    end
    return idx
end


function isInMaxComponent(C :: Array{Component,1}, u = nothing)
    if u != nothing
        idx = maxComponentIdx(C)
        if in(u, C[idx].nodemap) == true
            println("MAX VALUE $u IS IN COMPONENT $idx, WITH SIZE = ", C[idx].size)
        end
    end
end

function printDiagnostics(C :: Array{Component,1}, B :: Bridges)
    println("Bridges:")
    printBridges(B)
    println("Total components:", length(C))
    for (idx, c) in enumerate(C)
        print("$idx")
        printComponent(c)
    end
    println("End of printDiagnostics!\n")
end

function runlocalApprox(C :: Array{Component,1})
    for (idx, c) in enumerate(C)
        t = time()
        #print("Approximating component $idx ...")
        c.distances = localApprox(c, w)
        #println(" done")
        #println("calculate component time:", time() - t, "(s)")
    end
    return time()-t
end


function minDistanceInSet(nodes :: Array{Int64,1}, distances :: Array{Float64,1})
    if length(nodes) > length(distances)
        println("Error: array sizes do not match")
        exit()
    end
    min :: Float64 = Inf
    mu :: Int64 = 0
    for (idx,u) in enumerate(nodes)
        #logw(w," Distance of link $u = ", distances[idx])
        if distances[idx] < min
            min = distances[idx]
            mu = u
        end
    end
    return mu, min
end

function evaluateMinDistance(distances :: Array{Float64,1}, min :: Float64)
    total :: Int64 = 0 
    for i in distances
        if i < min
            total +=1
        end
    end
    return total
end




function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream, solut = nothing)
    start_time = time()
    u :: Int64 = 0
    n = A.n
    B = Bridges
    A, B = extractBridges(A)
    println("### extracting bridges time: ", time()- start_time, "(s)") 
    t = time()
    C = Array{Component,1}
    C = buildComponents(A, B)
    count = length(C)
    #println("-- Total components number:", count)
    println("### creating components time: ", time()- t, "(s)") 
    #isSpecialNode(B.edges, B.core2nodes, solut)
    #isInMaxComponent(C, solut)
    #printDiagnostics(C, B)
    return 
    if count == 1
        c = C[1]
        t = time()
        print("Approximating component 1 ... ")
        c.distances = localApprox(c, w)
        println(" done")
        println("calculate component time:", time() - t, "(s)")
        ### distances are basically kept on one dim array
        c.distances = sum(c.distances,2)
        cf = calculateCF(c.distances, n, c.nc)
        u = c.nodemap[indmax(cf)]
        logw(w,"node with argmax{c(", u , ")} = ", maximum(cf))
    else
        for (idx, c) in enumerate(C)
            t = time()
            print("Approximating component $idx ... ")
            c.distances = localApprox(c, w)
            println(" done")
            println("calculate component time:", time() - t, "(s)")
            lkdistances = sum(c.distances,1)[:,2:end]
            tldistances = sum(c.distances,2)[:]
            mu, min = minDistanceInSet(c.link, tldistances)
            rate = evaluateMinDistance(tldistances, min)
            println("for link: $mu, ", 100*rate/c.nc, "% nodes have smaller er!" ) 
            logw(w,"node with argmin{c(", c.nodemap[indmin(tldistances)], ")} = ", minimum(tldistances))
            ### TODO: can I calculate a link version calculation that saves time for JLT???
        end
        println("here!!")
        t = time()
        newcomp = B.comp
        newedges = B.edges
        n :: Int64 = length(newcomp)
        # following line: improves perf? TODO: check
        cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
        cA = contractAdjGraph(newedges, newcomp, C, count)
        #println(full(cA))
        # following 2 lines: improves perf? TODO: check
        dist2 = zeros(Float64,count,count)
        path2 = zeros(Int64,count,count)
        dist2, path2 = shortestContractPaths(cA, count, newedges, newcomp)
        #       println("path2:", path2)        
        #       println("dist2:", dist2)
        sizes = zeros(Int64,count)
        sizes = compRealSizes(C, count)
  #      println("Sizes:", sizes)
        distcomp = zeros(Float64,count)
        distcomp = compContractLocalDists(C, count, path2, dist2, sizes)
        #     println("Final distcomp:", distcomp)
        updateLocalDists(C, sizes,newedges, newcomp)
        fdistance, fnodes = aggregateLocalDists(C, distcomp)
        #fdistance, fnodes = aggregateDistances(C, count, path2, dist2, sizes, newedges, newcomp)
        #println(fdistance)
        #println(fnodes)
        cf = calculateCF(fdistance, A.n,length(fdistance))
        u = fnodes[indmax(cf)]
        logw(w,"\t node with argmax{c(", u , ")} = ", maximum(cf))
        println("updating distances ",time()- t, "(s)")
        println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
      end
    println("I am here!!")
    return u
end

### TODO: bug in following func
function wrcfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream, solution = nothing)
    println("****** Running (fast) approx ******")

    u  = cfcAccelerate(A, w, solution)
    if solution != nothing
        if solution != u
            logw(w,"\t THIS APPROX RESULT IS DIFFERENT THAN OTHERS (OR THE EXACT SOL = ", solution,")")
        end        
    end
end



function cfcAccelerate2(A:: SparseMatrixCSC{Float64}, w :: IOStream)
    start_time = time()
    u :: Int64 = 0
    B = Bridges
    A, B = extractBridges(A)
    println("### extracting bridges time: ", time()- start_time, "(s)")
    println("Bridges:")
    printBridges(B)
    t = time()
    C = Array{Component,1}
    
    C = buildComponents(A, B)
    count = length(C)
    println("### creating components time: ", time()- t, "(s)") 
    #printDiagnostics(C, B)
    println("Bridges:")
    printBridges(B)

    if count == 1
        c = C[1]
        t = time()
        print("Approximating component 1 ... ")
        c.distances = localApprox(c, w)
        println(" done")
        println("calculate component time:", time() - t, "(s)")
        ### distances are basically kept on one dim array
        c.distances = sum(c.distances,2)
        cf = calculateCF(c.distances, A.n, c.nc)
        u = c.nodemap[indmax(cf)]
        logw(w,"node with argmax{c(", u , ")} = ", maximum(cf))
    else
        solve_time = time()
        for (idx, c) in enumerate(C)
            t = time()
            print("Approximating component $idx ... ")
            c.distances = localApprox(c, w)
            println(" done")
            println("calculate component time:", time() - t, "(s)")
            lkdistances = sum(c.distances,1)[:,2:end]
            tldistances = sum(c.distances,2)[:]
            return
            #println("Distances: ",tldistances)
            #mu, min = minDistanceInSet(c.link, tldistances)
            #rate = evaluateMinDistance(tldistances, min)
            #println("for link: $mu, ", 100*rate/c.nc, "% nodes have smaller er!" ) 
            #logw(w,"node with argmin{c(", c.nodemap[indmin(tldistances)], ")} = ", minimum(tldistances))
            cf = calculateCF(tldistances, A.n, c.nc)
            u = c.nodemap[indmax(cf)]
            logw(w,"node with argmax{c(", u , ")} = ", maximum(cf))
            ### TODO: can I calculate a link version calculation that saves time for JLT???
        end
        println("### solving time: ", time()- solve_time, "(s)") 
        t = time()
        newcomp = B.comp
        newedges = B.edges 
        n :: Int64 = length(newcomp)
        # TODO: check if allocating memory for newly defined variables helps performances
        cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
        cA = contractAdjGraph(newedges, newcomp, C, count)
        #println(full(cA))
        return 
        dist2 = zeros(Float64,count,count)
        path2 = zeros(Int64,count,count)
        dist2, path2 = shortestContractPaths(cA, count, newedges, newcomp)
        #       println("path2:", path2)        
        #       println("dist2:", dist2)
        sizes = zeros(Int64,count)
        sizes = compRealSizes(C, count)
  #      println("Sizes:", sizes)
        distcomp = zeros(Float64,count)
        distcomp = compContractLocalDists(C, count, path2, dist2, sizes)
        #     println("Final distcomp:", distcomp)
        updateLocalDists(C, sizes,newedges, newcomp)
        fdistance, fnodes = aggregateLocalDists(C, distcomp)
        #fdistance, fnodes = aggregateDistances(C, count, path2, dist2, sizes, newedges, newcomp)
        #println(fdistance)
        #println(fnodes)
        cf = calculateCF(fdistance, A.n,length(fdistance))
        u = fnodes[indmax(cf)]
        logw(w,"\t node with argmax{c(", u , ")} = ", maximum(cf))
        println("updating distances ",time()- t, "(s)")
        println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
      end
    return u
end
