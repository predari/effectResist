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



#### beware of merges!
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
