include("structures.jl")
using Laplacians
using LightGraphs


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

    ## TODO: change to: for (idx, u) in enumerate(bdry)

    for (idx, u) in enumerate(bdry)
        #println("$idx $u")
        findin()
        cf[u,1] = sumer2 + n*er2[u] -2er[u]*sumer
        for i in 1:n
            cf[i,idx] = er2[i] + er2[u] -2er[i]*er[u]
        end
    end
    
    sumdeler2 = sum(delextnode(er2, bdry, bdryc))
    sumdeler = sum(delextnode(er, bdry, bdryc))
    for i in 1:n
        if (i in bdry) == false
            cf[i,1] = sumdeler2 + (n-1)*er2[i] -2er[i]*sumdeler
        end
    end
    return cf
end

#### TODO: check if I have matrix vector multiplications in here and perform it directly
function LinvDistance(c::Component ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    start_time = time()
    f = approxCholLap(c.A,tol=1e-5);

    nodes = c.nodemap
    #println(nodes)
    bdry = c.bdry
    external = c.external
    link = c.link
    sz = c.linkc

    n = c.nc
    k = round(Int, JLfac*log(n)) # number of dims for JL    
    U = wtedEdgeVertexMat(c.A)
    m = size(U,1)
    er = zeros(n)
    er2 = zeros(n)
    ## first position for Sigma(d(u,v: where v not in bdry)). For d(u,bdry)
    ## position in cf is the same as in the bdry array
    cf = zeros(n, sz + 1 )
    
    for i = 1:k # q 
        r = randn(m) 
        ur = U'*r 
        v = zeros(n)
        v = f(ur[:])
        er.+= v./k
        er2.+= v.^2/k
    end
    println("solving approxCholLap time:", time() - start_time, "(s)")
    
    t = time()
    sumer = sum(er)    
    sumer2 = sum(er2)


    #println("bdry:",bdry)
    l3c_idx = findin(nodes, link)
    l2c_idx = findin(nodes, bdry)
    #println("local:",l2c_idx)
    #println("external:",external)
    
    #println("l2c_idx (local numb):", l2c_idx)
    #println("l1c_idx (local numb):", l1c_idx)

    for (idx, u) in enumerate(l3c_idx)
        cf[:,idx + 1] .= er2 + er2[u] .- 2er*er[u]
        # for i in 1:n
        #     cf[i, idx + 1] = er2[i] + er2[u] -2er[i]*er[u]
        # end
        sumer += - er[u]
        sumer2 += - er2[u]
    end

    #multiplier = external[findin(bdry,getindex(nodes,l1c_idx))]
    #println("l1c_idx (local numb):", l1c_idx)
    #println("multiplier :", multiplier)

    cf[:,1] .= (n-sz)*er2 + sumer2 .-2er*sumer
    # for u in 1:n
    #     cf[u,1] =  sumer2 + (n-sz)*er2[u] -2er[u]*sumer
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
        cf[:,1] = cf[:,1] .+ (er2 + er2[v] .-2er*er[v]) * external[idx] + external[idx]
    end

    println("updating distances time:", time() - t, "(s)")
    
    return cf
end




function LinvDistance(a::SparseMatrixCSC{Float64}, bdry::Array{Int64,1}, bsize:: Int64, nodes::Array{Int64,1} ; ep=0.3, matrixConcConst=4.0, JLfac=200.0)
    if 0 in bdry
        println("Error: node index cannot be zero!");
        exit(0)
    end
    f = approxCholLap(a,tol=1e-5);
    n = size(a,1)
    k = round(Int, JLfac*log(n)) # number of dims for JL    
    U = wtedEdgeVertexMat(a)
    m = size(U,1)
    println(bdry)
    println(n)
    er = zeros(n)
    er2 = zeros(n)
    cf = zeros(n, bsize + 1)
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

    println("old er2", er2)
    sumer2 = sum(er2)
    sumer = sum(er)
    println("size er2: ",size(er2,1))
    ## TODO: change to: for (idx, u) in enumerate(bdry)
    lu_array = Array{Int64,1}(bsize)
    for (idx, u) in enumerate(bdry)
        tmp = findin(nodes, u)
        lu = tmp[1] ## TODO: why I need to convert
        println("idx:$idx global-u:$u")
        lu_array[idx] = lu
        #cf[lu, 1] = sumer2 + n*er2[lu] -2er[lu]*sumer
        println("cf of local-link $lu:",cf[lu,1])
        for i in 1:n
            cf[i, idx + 1] = er2[i] + er2[lu] -2er[i]*er[lu]
        end
        ## or remove node u from er2 and er here already
        ## no need to store in this case.
    end
    println("Bdry nodes (local numb):", lu_array)
    # sumdeler2 = sum(er2) #delextnode(er2, lu_array, bsize)
    # sumdeler = sum(er2) #delextnode(er, lu_array, bsize)

    for i in lu_array
        sumer = sumer -er[i]
        sumer2 = sumer2 -er2[i]
    end
#    println("new er2", ll)
#    println("old er2", er2)
#    sumdeler2 = sum(ll)
#    sumdeler = sum(delextnode(er, lu_array, bsize))
#    if (n - bsize) != size(ll,1)
#        println(n - bsize, size(ll,1))
#        println("ERROR: something is wrong!!")
#    end
    for i in 1:n
        #if (i in lu_array) == false
            # (n-1) before
            cf[i,1] = sumer2 + (n-bsize)*er2[i] -2er[i]*sumer
     end
    return cf
end

