include("logw.jl")
include("graph.jl")
include("lapl.jl")
include("Linvdiag.jl")
include("Lpartinv.jl")
include("Lpartinv2.jl")
include("Lpartinv3.jl")
include("appxInvTrace.jl")

using Laplacians

function delnode2(L:: SparseMatrixCSC{Float64}, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function approx(A, L, k, stretch, flag, w = open("log.txt", "w")) #v- the vertex chosed; G- graph; k- the number of edge added

    n = A.n
    ans = zeros(2,1)
    t = time()
    #rank = zeros(Int64,k,1)
    #rank = zeros(Int64,k,n)
    rank = zeros(Float64,k,n)
    L
    #logw(w,"\t k=1 starts")
    start_time=time()

    t=time()
    Ldele = LinvdiagSS(A;JLfac=20)
   # logw(w,"Find 1st u with argmax{c(u)}")
   # logw(w,"Ldele: ")
   # print(Ldele)
    # logw(w,"LaplSolve total time: ",time()-t, " (s)")
    
    u = indmin(Ldele)
    #rank[1] = u
    sl = sortperm(Ldele)
    logw(w,"\t node ranking (min -> max):", sl[1:10])
    sv = sort(Ldele)
    logw(w,"\t cfcc ranking (min -> max):", sv[1:10])
    rank[1,:] = sl
    #logw(w,"")
    
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)

    elapsed_time = time()-start_time

    start_time = time()
    for dep = 1:k-1
        #logw(w,"\t k=",dep+1," starts")

        t=time()
        f = approxCholSddm(L,tol=1e-5, verbose=false);
        #logw(w,"\t SDDM preconditioning time: ",time()-t, " (s)")
        t=time()
        Ldele = LpartinvSS(n - dep, f;JLfac=20)

        #logw(w,"\t SDDM1 solve time: ",time()-t, " (s)")
        t=time()
        Ldele2 = Lpartinv2SS(n - dep, f, A;JLfac=20)
        #logw(w,"\t SDDM2 solve time: ",time()-t, " (s)")
        t=time()
        Ldele3 = Lpartinv3SS(L, f;JLfac=20)
        #logw(w,"\t SDDM3 solve time: ",time()-t, " (s)")

        Ldele ./= (Ldele2 .+ Ldele3)


        u = indmax(Ldele)
        sl = sortperm(Ldele[:,1], rev=true)
        logw(w,"\t node ranking (min -> max):", sl[1:10])
        sv = sort(Ldele[:,1], rev=true)
        logw(w,"\t cfcc ranking (min -> max):", sv[1:10])
        #rank[dep+1] = u
        x = zeros(Int64,dep)
        append!(sl,x)
        rank[dep+1,:] = sl
        
        L = delnode2(L,u,n-dep)
        A = delnode2(A,u,n-dep)
    end

    elapsed_time += time()-start_time
    logw(w,"\t approx_elapsed_time = ",elapsed_time, "s")
    ans[1] = elapsed_time
    #logw(w,"\t calculating CFCC of the solution returned by approx greedy..")
    # originally (n > 30000)
    # the following calculation takes time. Avoid it with flag = 0
    if(flag == 1)    
        ans[2] = (n > 20000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    end
    if (stretch != 0)
        ans[2] = stretch * ans[2]
    end
    logw(w,"\t CFCC achieved by approx greedy: ", ans[2])
    return ans, rank
end
#######
# WARNING: do not use it! Takes a lot of time, and returns not valuable info.
#######
# returns the corresponding approximation cf values for a given set of nodes.
# I cannot really compare the augtree method because of the stretch, disorting the values
# (although the ranking of the nodes can be used to compare the methods).
function check_approx_values(A, k, nodes::Array{Int64,2}, w = open("log.txt", "w"))

    n = A.n
    ans = zeros(3,1)
    t = time()
    rank = zeros(Int64,k,1)
    cf = zeros(k,1)
    L = lap(A)
    start_time=time()

    t=time()
    Ldele = LinvdiagSS(A;JLfac=20)
    u = nodes[1]
    cf[1] = Ldele[u]
    
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)

    elapsed_time = time()-start_time

    start_time = time()
    for dep = 1:k-1

        t=time()
        f = approxCholSddm(L,tol=1e-5, verbose=false);
        t=time()
        Ldele = LpartinvSS(n - dep, f;JLfac=20)
        t=time()
        Ldele2 = Lpartinv2SS(n - dep, f, A;JLfac=20)
        t=time()
        Ldele3 = Lpartinv3SS(L, f;JLfac=20)
        Ldele ./= (Ldele2 .+ Ldele3)

        u = nodes[dep + 1]
        cf[dep + 1] = Ldele[u,1]
        rank[dep+1] = u

        L = delnode2(L,u,n-dep)
        A = delnode2(A,u,n-dep)
    end
    println("cf:", cf)
    return (n > 20000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
end
