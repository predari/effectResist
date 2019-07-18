
include("logw.jl")
include("LinvDistance.jl")
include("ShortestPaths.jl")
include("structures.jl")
include("bridge.jl")
include("components.jl")
include("compDistances.jl")

include("sampling.jl")
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
    logw(w,"****** Running fast approx ******")
    distances = LinvDistance(c)
    return distances
end


function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream)
    start_time = time()
    n = A.n
    #computeCore2Bridges(A)
    B = Bridges
    A, B = extractBridges(A)
    t = time()
    C = Array{Component,1}
    C = buildComponents(A, B)
    count = length(C)
    for c in C
        print(c.size,"-",length(c.nodemap)," ")
    end
    println("creating components time: ", time()- t, "(s)") 
    t = time()
    for (idx, c) in enumerate(C)
        print("Approxing component $idx ...")
        c.distances = localApprox(c, w)
        println(" done")
    end
    # println("Bridges:")
    # printBridges(B)
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
        #println("components:", B.comp)
        #println("edges:", B.edges)
        #println("size:",length(B.comp))
        ## strip core1nodes in order to remove internal paths
        ##### stripcore1nodes is in components.jl
        newcomp, newedges = stripcore1nodes(B.comp,B.edges)
        #println("After striping nodes1!")
        #println("components:", newcomp)
        #println("size:",length(newcomp))
        
        #println("edges:", newedges)
        n :: Int64 = length(newcomp)
        # following line: improves perf? TODO: check
        cA :: SparseMatrixCSC{Float64} = spzeros(n,n)
        cA = contractAdjGraph(newedges, newcomp, C, count)
        #println(cA)
        # following 2 lines: improves perf? TODO: check
        dist2 = zeros(Float64,count,count)
        path2 = zeros(Int64,count,count)
        dist2, path2 = shortestContractPaths(cA, count, newedges, newcomp)
        #println("path2:", path2)        
        #println("dist2:", dist2)
        sizes = zeros(Int64,count)
        sizes = compRealSizes(C, count)
        #println("Sizes:", sizes)
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



function cfcAccelerate(A:: SparseMatrixCSC{Float64}, w :: IOStream, maxcf :: Int64)
    start_time = time()
    n = A.n
    B = Bridges
    A, B = extractBridges(A)
    t = time()
    C = Array{Component,1}
    C = buildComponents(A, B)
    count = length(C)
    println("number of components is = ", count)
    # println("Bridges:")
    # printBridges(B)
    # println("Components: $count")
    # for (idx, c) in enumerate(C)
    #     print("$idx")
    #     printComponent(c)
    # end
    #println("maxcf = $maxcf")
    # for c in C
    #     print(c.size,"-",length(c.nodemap)," ")
    # end
    # println()
    for node in B.edges
        if maxcf == node
            println("MAX VALUE $node IS PART OF LINK NODES!!!")
            break;
        end
    end
    for node in B.core2nodes
        if maxcf == node
            println("MAX VALUE $node IS PART OF CORE2NODES NODES!!!")
            break;
        end
    end
    for (idx, c) in enumerate(C)
        if findin(c.nodemap, maxcf) != Array{Int64,1}(0)
            println("max is in component $idx with length=", length(c.nodemap))
            break;
        end
    end
    println("creating components time: ", time()- t, "(s)")

    # for (idx, c) in enumerate(C)
    #     t = time()
    #     #print("Approxing component $idx ...")
    #     c.distances = localApprox(c, w)
    #     println("calculate component $idx (", c.nc, ") time:", time() - t, "(s)")
    #     #println(" done")
    # end
    # println(" solving core1 time : ", time()- t, "(s)")
    if count == 1
        c = C[1]
        t = time()
        c.distances = localApprox(c, w)
        println("calculate component 1 (", c.nc, ") time:", time() - t, "(s)")
        c.distances = sum(c.distances,2)
        cf = calculateCF(c.distances, n, c.nc)
        logw(w,"\t node with argmax{c(", c.nodemap[indmax(cf)], ")} = ", maximum(cf))
    else
        maxcmpsize = 0
        maxcmpid = 0
        for (idx, c) in enumerate(C)
            if c.size > maxcmpsize
                maxcmpsize = c.size
                maxcmpsize = idx
            end
            t = time()
            c.distances = localApprox(c, w)
            println("calculate component $idx (", c.nc, ") time:", time() - t, "(s)")
            linkdistance = sum(c.distances,1)[:,2:end]
            localdistances = sum(c.distances,2)
            #cf = calculateCF(c.distances, n, c.nc)
            logw(w,"\t Locally: node with argmin{c(", c.nodemap[indmin(localdistances)], ")} = ", minimum(localdistances))
            #println(c.distances)
            for (idx,u) in enumerate(c.link)
                #logw(w,"\t Link $u = ", c.distances[findin(c.nodemap, u)])
                tmp = 0
                for i in c.distances
                    if i < getindex(localdistances[findin(c.nodemap, u)])
                        tmp +=1
                    end
                end
                #println("For link $u ", 100*tmp/c.nc, "% nodes have smaller effective resistance!" )
            end
            pivots = Array{Int64,1}
            pv :: Int64 = 10
            if c.nc < pv
                println("WARNING: number of pivots is the same as the subgraph.")
                pv = c.nc
            end
            if c.nc > 1
                f = approxCholLap(c.A,tol=1e-5);
            else
                f = nothing
            end
            pivots =samplePivots(c.nc, pv)
            t = time()
            SamplingDistAll(c, pivots, f)
            println("calculate sampling time (all):", time() - t, "(s)")
            t = time()
            distances = SamplingDistLink(c, pivots,f)
            println("calculate sampling time (links):", time() - t, "(s)")
            tmp = 0
            for (idx,u) in enumerate(c.link)
                #logw(w,"\t Link ", u," = ", distances[idx])
                for i in c.distances
                    if i < getindex(distances[idx])
                        tmp +=1
                    end
                end
                #println("For link $u ", 100*tmp/c.nc, "% nodes have smaller effective resistance!" )
                    
            end
        end

        t = time()
        #println("components:", B.comp)
        #println("edges:", B.edges)
        ## strip core1nodes in order to remove internal paths
        #newcomp, newedges = stripcore1nodes(B.comp,B.edges)
        #println("After striping nodes1!")
        newcomp =B.comp
        newedges =B.edges
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
        logw(w,"\t node with argmax{c(", fnodes[indmax(cf)], ")} = ", maximum(cf))
        println("updating distances ",time()- t, "(s)")
        println("TOTAL CFC TIME IS: ", time() - start_time, "(s)")
    end
end

