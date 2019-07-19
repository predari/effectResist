include("structures.jl")
using StatsBase ## for pivots
using Laplacians

function samplePivots(n :: Int64, pv :: Int64  )
    pivots = Array{Int64,1}
    if  n < pv
        println("WARNING: number of pivots is smaller than number of nodes!")
        exit()
    end
    pivots = StatsBase.sample(1:n, pv, replace = false)
    return pivots
end



#### beware of merges!
function SamplingDistLink(c :: Component, pivots :: Array{Int64,1}, lapSolver)
    start_time = time()
    n = c.nc
    sz = c.linkc
    if n == 1
        return zeros(sz)
    end
    nodes = c.nodemap
    external = c.external
    link = c.link


    l3c_idx = findin(nodes, link)
    pv = length(pivots)
    cf = zeros(n)
    # println("nodes:",nodes)
    # println("link:",link)
    # println("l3c_idx:",l3c_idx)
    # println("pivots:",pivots)
    
    b = zeros(n)
    er = zeros(sz)
    for (idx1, u) in enumerate(l3c_idx)
        b[u] = 1.0
        for (idx2, s) in enumerate(pivots)
            if b[s] == 0
                b[s] = -1.0
            end
            v = zeros(n)
            v = lapSolver(b)
            er[idx1] += v[u] - v[s]
            b[s] = 0.0
        end
    end
    
    return er
end






function SamplingDistAll(A::SparseMatrixCSC{Float64}, pivots :: Array{Int64,1})
    start_time = time()
    f = approxCholLap(A,tol=1e-5);
    n = A.n
    if n == 1
        return zeros(n)
    end
    pv = length(pivots)
    cf = zeros(n)
    for u in 1:n
        b = zeros(n)
        b[u] = 1.0
        for s in enumerate(pivots)
            if b[s] == 0
                b[s] = -1.0
            end
            v = zeros(n)
            v = f(b)
            cf[u] += v[u] - v[s]
            b[s] = 0.0
        end
    end   
    return cf
end






function SamplingDistAll(c :: Component, pivots :: Array{Int64,1}, lapSolver)
    start_time = time()
    n = c.nc
        if n == 1
        return zeros(n)
    end
    nodes = c.nodemap
    pv = length(pivots)
    cf = zeros(n)
    for (idx1, u) in enumerate(nodes)
        b = zeros(n)
        b[idx1] = 1.0
        for (idx2, s) in enumerate(pivots)
            if b[s] == 0
                b[s] = -1.0
            end
            v = zeros(n)
            v = lapSolver(b)
            cf[idx1] += v[idx1] - v[s]
            b[s] = 0.0
        end
    end   
    return cf
end



#### beware of merges!
function SamplingDistLink(c :: Component, pv :: Int64 ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    start_time = time()
    sz = c.linkc
    n = c.nc
    if n == 1
        return zeros(pv, sz)
    end    
    nodes = c.nodemap
    bdry = c.bdry
    external = c.external
    link = c.link

    l3c_idx = findin(nodes, link)
    l2c_idx = findin(nodes, bdry)
    
    #println("size of link = ", sz)
    #println("local link num= ", l3c_idx)
    
    t = time()
    f = approxCholLap(c.A,tol=1e-5);
    println(" f = CholLap time: ", time()- t, "(s)")
    if  n < pv
        println("WARNING: number of pivots is smaller than number of nodes!")
        exit()
    end
    pivots = StatsBase.sample(1:n, pv, replace = false)
    cf = zeros(pv, sz)

    b = zeros(n)
    for (idx1, u) in enumerate(l3c_idx)
        b[u] = 1.0
        for (idx2, s) in enumerate(pivots)  
            b[s] = -1.0            
            v = zeros(n)
            v = f(b)
            cf[idx2,idx1] = v[u] - v[s]
            b[s] = 0.0
        end
    end   
    return cf
end

#### beware of merges!
function SamplingDistAll(c :: Component, pv :: Int64 ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    start_time = time()
    sz = c.linkc
    n = c.nc
    if n == 1
        return zeros(n)
    end
    nodes = c.nodemap
    t = time()
    f = approxCholLap(c.A,tol=1e-5);
    println(" f = CholLap time: ", time()- t, "(s)")
    if  n < pv
        println("WARNING: number of pivots is smaller than number of nodes!")
        exit()
    end
    pivots = StatsBase.sample(1:n, pv, replace = false)
    cf = zeros(n)

    b = zeros(n)
    for (idx1, u) in enumerate(nodes)
        b[idx1] = 1.0
        for s in pivots  
            b[s] = -1.0            
            v = zeros(n)
            v = f(b)
            cf[idx1] += v[idx] - v[s]
            b[s] = 0.0
        end
    end   
    return cf
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
    println("-- Total components number:", count)
    println("### creating components time: ", time()- t, "(s)") 
    isSpecialNode(B.edges, B.core2nodes, solut)
    isInMaxComponent(C, solut)
    printDiagnostics(C, B)

    #solveTime = runlocalApprox(C)
    # println(" solving time : ", solvetime, "(s)")
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
    end
end

function samplingCompApprox(C :: Component, pv:: Int64, solution = nothing)
    start_time = time()
    u :: Int64 = 0
    pivots = Array{Int64,1}
    pv :: Int64 = 10
    for (idx, c) in enumerate(C)
        t = time()
        if c.nc < pv
            println("WARNING: setting the number of pivots to be same as the subgraph.")
            pv = c.nc
        end
        if c.nc > 1
            f = approxCholLap(c.A,tol=1e-5);
        else
            f = nothing
        end
        pivots =samplePivots(c.nc, pv)
        println(" setting up pivots time : ", time()- t, "(s)")
        t = time()
        c.distances = SamplingDistAll(c, pivots, f)
        println("calculate sampling time (all):", time() - t, "(s)")
        mu, min = minDistanceInSet(c.link, c.distances)
        rate = evaluateMinDistance(c.distances, min)
        println("for link: $mu, ", 100*rate/c.nc, "% nodes have smaller er!" ) 
        t = time()
        distances = SamplingDistLink(c, pivots,f)
        println("calculate sampling time (links):", time() - t, "(s)")
    end
end


function samplingApprox(A:: SparseMatrixCSC{Float64}, w :: IOStream, pv:: Int64, solution = nothing)
    start_time = time()
    u :: Int64 = 0
    n = A.n
    pivots = Array{Int64,1}
    pivots =samplePivots(n, pv)
    println(" setting up pivots time : ", time()- start_time, "(s)")
    t = time()
    distances = SamplingDistAll(A, pivots)
    println("calculate sampling time:", time() - t, "(s)")
    cf = calculateCF(distances, n, n)
    u = indmax(cf)
    logw(w,"node with argmax{c(", u , ")} = ", maximum(cf))
end


function wrsamplingApprox(A:: SparseMatrixCSC{Float64}, w :: IOStream, pv:: Int64, solution = nothing)
    u  = cfcAccelerate(A, w, solution)
    if solution != nothing
        if solution != u
            logw(w,"\t THIS APPROX RESULT IS DIFFERENT THAN OTHERS (OR THE EXACT SOL = ", solution,")")
        end        
    end
end
