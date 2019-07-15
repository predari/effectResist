include("structures.jl")
using Laplacians
using LightGraphs
using StatsBase

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

#### beware of merges!
function SamplingDistLink(c :: Component, pivots :: Array{Int64,1}, lapSolver)
    start_time = time()
    
    nodes = c.nodemap
    external = c.external
    link = c.link
    sz = c.linkc
    n = c.nc
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
    nodes = c.nodemap
    n = c.nc
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
    
    nodes = c.nodemap
    bdry = c.bdry
    external = c.external
    link = c.link
    sz = c.linkc
    n = c.nc

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
    
    nodes = c.nodemap
    sz = c.linkc
    n = c.nc

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



#### beware of merges!
function LinksLinvDist(c::Component ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    start_time = time()
    
    nodes = c.nodemap
    bdry = c.bdry
    external = c.external
    link = c.link
    sz = c.linkc
    n = c.nc
    println(nodes)
    println(link)

    l3c_idx = findin(nodes, link)
    l2c_idx = findin(nodes, bdry)
    #println("local link numbering ", l3c_idx)

    println("size of link = ", sz)
    println("link = ", l3c_idx)
    t = time()
    f = approxCholLap(c.A[[l3c_idx],[l3c_idx]],tol=1e-5);
    f = cholLap(a,tol=1e-5);
    println(" f = CholLap time: ", time()- t, "(s)")
    k = round(Int, JLfac*log(n))     
    U = wtedEdgeVertexMat3(c.A,l3c_idx)
    m = size(U,1)
    l = size(U,2)
    println("m = ", m, " nnz(A) = ", nnz(c.A)/2)
    println("l = ", l, " sz = ", sz)
    er = zeros(l)
    er2 = zeros(l)
    cf = zeros(l, sz + 1)
 
    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(l)
        v = f(ur[:])
        er.+= v./k
        er2.+= v.^2/k
    end
    
    for (idx, u) in enumerate(l3c_idx)
        for i in 1:n
            cf[i, idx + 1] = er2[i] + er2[u] -2er[i]*er[u]
        end
    end
   
    println("calculating only for link time (without creating of matrix):", time() - t, "(s)")
    return cf
end


### check for correctness!!!
function wtedEdgeVertexMat3(mat::SparseMatrixCSC, nodes ::Array{Int64,1})
    ## return U (m,n) only for link nodes --> U (m,l)    
    (ai,aj,av) = findnz(triu(mat,1))

    m = length(ai)
    n = size(mat)[2] 
    v = av.^(1/2)
    return (sparse(collect(1:m),ai,v,m,n) - sparse(collect(1:m),aj,v,m,n))[:,nodes]
end



function LinvDistance(A::SparseMatrixCSC{Float64}, bdry::Array{Int64,1},  external :: Array{Array{Int, 1}, 1}, nodes::Array{Int64,1} ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    start_time = time()
    f = approxCholLap(A,tol=1e-5);
    n = A.n
    k = round(Int, JLfac*log(n)) # number of dims for JL    
    U = wtedEdgeVertexMat(A)
    m = size(U,1)
    println("m = ", m, " nnz(A) = ", nnz(A)/2)
    er = zeros(n)
    er2 = zeros(n)
    cf = zeros(n)
    
    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(n)
        v = f(ur[:])
        er.+= v./k
        er2.+= v.^2/k
    end

    t = time()
    sumer = sum(er)    
    sumer2 = sum(er2)
    cf[:,1] .= (n)*er2 + sumer2 .-2er*sumer
    
    l2c_idx = findin(nodes, bdry)
    for (idx, v) in enumerate(l2c_idx)
        m1 :: Int64 = 0
        m2 :: Int64 = 0
        for (i, ex) in enumerate(external[idx])
            m1 += ex 
            m2 += i*ex
        end
        cf[:,1] = cf[:,1] .+ (er2 + er2[v] .-2er*er[v]) * m1 + m2
    end
    return cf
end

#### TODO: check if I have matrix vector multiplications in here and perform it directly
function LinvDistance(c::Component ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    start_time = time()
    f = approxCholLap(c.A,tol=1e-5);

    nodes = c.nodemap
    bdry = c.bdry
    external = c.external
    link = c.link
    sz = c.linkc

    n = c.nc
    k = round(Int, JLfac*log(n)) # number of dims for JL    
    U = wtedEdgeVertexMat(c.A)
    m = size(U,1)
    println("m = ", m, " nnz(A) = ", nnz(c.A)/2)
    er = zeros(n)
    er2 = zeros(n)
    ## first position for Sigma(d(u,v: where v not in bdry)). For d(u,bdry)
    ## position in cf is the same as in the bdry array
    cf = zeros(n, sz + 1)
    
    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(n)
        v = f(ur[:])
        er.+= v./k
        er2.+= v.^2/k
    end
    #println("solving approxCholLap time:", time() - start_time, "(s)")
    #er = sqrt.(er2)
    # println("er",er)
    # println("er2",er2)
    # println("er2[2]",er2[2], "er[2] ",er[2])
    t = time()
    sumer = sum(er)    
    sumer2 = sum(er2)

    l3c_idx = findin(nodes, link)
    l2c_idx = findin(nodes, bdry)
#    println("l3c_idx (local link numb):", l3c_idx)
#    println("l2c_idx (local bdry numb):", l2c_idx)

    for (idx, u) in enumerate(l3c_idx)
        #cf[:,idx + 1] .= er2 + er2[u] .- 2er*er[u]
        for i in 1:n
            cf[i, idx + 1] = er2[i] + er2[u] -2er[i]*er[u]
        end
        sumer += - er[u]
        sumer2 += - er2[u]
    end
    #println(cf)

    #multiplier = external[findin(bdry,getindex(nodes,l1c_idx))]
    #println("l1c_idx (local numb):", l1c_idx)
    #println("multiplier :", multiplier)

    cf[:,1] .= (n-sz)*er2 + sumer2 .-2er*sumer
    # for u in 1:n
    #     #cf[u,1] =  sumer2 + (n-sz)*er2[u] -2er[u]*sumer
    #     ### I do the following for one extra time for u == v.
    #     ### To avoid do for v in setdiff(l1c_idx,u), but then
    #     ### the idx should be recomputed with findin etc.
    #     for (idx, v) in enumerate(l2c_idx) 
    #         #cf[u,1] = cf[u,1] + (er2[u] + er2[v] -2er[u]*er[v]) * multiplier[idx] + multiplier[idx]
    #         cf[u,1] = cf[u,1] + (er2[u] + er2[v] -2er[u]*er[v]) * external[idx] + external[idx]
    #     end
    # end
    #### Here, I calculate this value ss:(er2 + er2[v] .-2er*er[v]) for each v again!!
    #### I do that a first time in L144. I can avoid it by doing for those v:
    #### ss*(external[idx]+1) + external[idx] in L144
    for (idx, v) in enumerate(l2c_idx)
        m1 :: Int64 = 0
        m2 :: Int64 = 0
        for (i, ex) in enumerate(external[idx])
            m1 += ex 
            m2 += i*ex
        end
        #cf[:,1] = cf[:,1] .+ (er2 + er2[v] .-2er*er[v]) * external[idx] + external[idx]
        cf[:,1] = cf[:,1] .+ (er2 + er2[v] .-2er*er[v]) * m1 + m2
    end

    #println("updating distances time:", time() - t, "(s)")
    return cf
end


#### TODO: check if I have matrix vector multiplications in here and perform it directly
function LinvDistanceLinks(c::Component ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    start_time = time()
    
    nodes = c.nodemap
    bdry = c.bdry
    external = c.external
    link = c.link
    sz = c.linkc
    n = c.nc
    # rows = rowvals(c.A)
    # vals = nonzeros(c.A)
    # for u in 1:n
    #     for v in nzrange(c.A, u)
    #         print(vals[v]," ")
    #         c.A[u,rows[v]]
    #     end
    #     println()
    # end
    # println(full(c.A))
    # println("matrix:")
    # rows = rowvals(A)
    # vals = nonzeros(A)
    # for i = 1:k
    #     for j in nzrange(A, i)
    #         row = rows[j]
    #         val = vals[j]
    #         print(val," ")
    #     end
    #     println()
    # end
    # println(size(A,1), " ", size(A,2))

    l3c_idx = findin(nodes, link)
    l2c_idx = findin(nodes, bdry)
    #println("local link numbering ", l3c_idx)
    # vec = sparse(c.A[l3c_idx,:])
    # eye = zeros(Int64, size(vec,2),size(vec,1))
    # #println(eye)
    # for (idx,u) in enumerate(l3c_idx)
    #     eye[u,idx] = 1
    # end
    
    # #println("eye",eye)
    # A :: SparseMatrixCSC{Float64} = eye * vec
    # for i in 1:size(A,1)
    #     for j in 1:size(A,2)
    #         if A[i,j] == 1
    #             A[j,i] = 1
    #         end
    #     end
    # end
    t = time()
    println("size of link",sz)
    f = approxCholLap(c.A,tol=1e-5);
    println(" f = approxCholLap time: ", time()- t, "(s)")
    ### TODO: log n or smthing else?
    k = round(Int, JLfac*log(n))     
    U = wtedEdgeVertexMat2(c.A[l3c_idx,:])
    m = size(U,1)
    # if m is smaller than nnz(c.A)/2, then that's good!
    println("m = ", m, " nnz(A) = ", nnz(c.A)/2)
    er = zeros(n)
    er2 = zeros(n)
    cf = zeros(n, sz + 1)
 
    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(n)
        v = f(ur[:])
        #println(size(v,1), " ", size(v,2))
        er.+= v./k
        er2.+= v.^2/k
    end
    #println("solving approxCholLap time:", time() - start_time, "(s)")
    sumer = sum(er)    
    sumer2 = sum(er2)

    
    for (idx, u) in enumerate(l3c_idx)
        for i in 1:n
            cf[i, idx + 1] = er2[i] + er2[u] -2er[i]*er[u]
        end
        ## TODO: delete following two lines
        sumer += - er[u]
        sumer2 += - er2[u]
    end
    
    # for (idx, v) in enumerate(l2c_idx)
    #     m1 :: Int64 = 0
    #     m2 :: Int64 = 0
    #     for (i, ex) in enumerate(external[idx])
    #         m1 += ex 
    #         m2 += i*ex
    #     end
    #      cf[:,1] = cf[:,1] .+ (er2 + er2[v] .-2er*er[v]) * m1 + m2
    # end
    println("calculating only for link time (without creating of matrix):", time() - t, "(s)")
    return cf
end


#### TODO: check if I have matrix vector multiplications in here and perform it directly
function LinvDistanceLinks2(c::Component ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    start_time = time()
    nodes = c.nodemap
    bdry = c.bdry
    external = c.external
    link = c.link
    sz = c.linkc
    println("size of link = ",sz)
    n = c.nc

    l3c_idx = findin(nodes, link)
    l2c_idx = findin(nodes, bdry)

    t = time()
    f = approxCholLap(c.A,tol=1e-5);
    println(" f = approxCholLap time: ", time()- t, "(s)")
    ### TODO: log n or smthing else?
    k = round(Int, JLfac*log(n))     
    U = wtedEdgeVertexMat(c.A)
    m = size(U,1)
    println("m = ", m, " nnz(A) = ", nnz(c.A)/2, " n = ",size(U,1))
    U = U[:,l3c_idx]
    println("reduced size: m = ",size(U,1), " n = ", size(U,2))
    er = zeros(sz)
    
    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(sz)
        if size(ur,1) == 1
            v = ur
        else
            v = f(ur[:])
        end
        #println(size(v,1), " ", size(v,2))
        er.+= v.^2/k
    end    
    println("calculating only for link time (without creating of matrix):", time() - t, "(s)")
    return er
end


### check for correctness!!!
function wtedEdgeVertexMat2(mat::SparseMatrixCSC)
    ## this is only for collection of rows. not for symmetric matrices
    ##(ai,aj,av) = findnz(triu(mat,1))
    
    r = size(mat)[1]
    if r != 1       
        (ai,aj,av) = findnz(triu(mat,1))
    else
        (ai,aj,av) = findnz(mat)
    end
    m = length(ai)
    n = size(mat)[2] 
    v = av.^(1/2)
    return sparse(collect(1:m),ai,v,m,n) - sparse(collect(1:m),aj,v,m,n)
end

