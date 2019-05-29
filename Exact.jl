include("graph.jl")
include("lapl.jl")
include("logw.jl")


function delnode2(L, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function delnode(Linv, v, t)
    Linv2 = (Linv - Linv[:,v]*Linv[:,v]'/(Linv[v,v]))
    return Linv2[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function exact(G, k,  stretch, w :: IOStream) #G- graph; k- the number of vertices in set S

    n = G.n;
    ans = zeros(2,1);
    #rank = zeros(Int64,k,1)
    #rank = zeros(Int64,k,n)
    rank = zeros(Float64,k,n)
    L = lapl(G)

    start_time = time()
    Linv = mppinv(L)

    min = Linv[1,1]
    u = 1
    for i = 2:n
        if Linv[i,i] < min
            u = i
            min = Linv[i,i]
        end
    end
    #logw(w, "k = 1: \t Linv.n: ", size(Linv,1), "\t Linv.m: ", size(Linv,2))
    dL = diag(Linv,0)
    #print(dL)
    sl = sortperm(dL)
    logw(w,"\t node ranking (min -> max):", sl[1:10])
    sv = sort(dL)
    logw(w,"\t cfcc ranking (min -> max):", sv[1:10])
 
    #rank[1] = u
    rank[1,:] = sl

    for i = 1:n
        L2 = delnode2(L,i,n)
        Linv = inv(L2)
        println("L2",L2)
        println("Linv after deling node:",i)
        println(Linv)
        println("trace:", trace(Linv))
    end
    
    L2 = delnode2(L,u,n)
    Linv = inv(L2)

    elapsed_time = time()-start_time

    #logw(w,"k = 1 total time: ",time() - elapsed_time, " (s)")

    start_time = time()
    for dep = 1:k-1
        t = time()
        dec = [(sum(Linv[:,i].^2)/Linv[i,i]) for i = 1:n-dep]
        u = indmax(dec)
        # i need u with max cc
        sl = sortperm(dec, rev=true)
        logw(w,"\t node ranking (min -> max):", sl[1:10])
        sv = sort(dec, rev=true)
        logw(w,"\t cfcc ranking (min -> max):", sv[1:10])
        #rank[dep+1] = u
        x = zeros(Int64,dep)
        append!(sl,x)
        rank[dep+1,:] = sl
        
        Linv = delnode(Linv,u,n-dep)
        #logw(w,"k = ", dep + 1 ," total time: ", time() - t, " (s)")
    end

    elapsed_time+=time()-start_time
    logw(w,"\t exact_elapsed_time = ", elapsed_time, " (s)")
    ans[1] = elapsed_time
    #logw(w,"\t calculating CFCC of the solution returned by exact greedy..")
    #logw(w, "k = ", k, "\t Linv.n: ", size(Linv,1), "\t Linv.m: ", size(Linv,2))
    #println("diag Linv: ", diag(Linv,0))
    ans[2] = n/trace(Linv)
    if (stretch != 0)
        ans[2] = stretch * ans[2]
    end
    logw(w,"\t CFCC achieved by exact greedy: ", ans[2])
    #logw(w,"")
    return ans, rank
end


function check_exact_values(G, nodes::Array{Int64,1}, k, w :: IOStream) #G- graph

    n = G.n;
    L = lapl(G)
    l = size(nodes,1)
    if l!=k
        exit()
    end
    cf = zeros(k,1)
    Linv = mppinv(L)    
    dL = diag(Linv,0)
    u = nodes[1]
    cf[1] = dL[u]
    L2 = delnode2(L,u,n)
    Linv = inv(L2)

    for dep = 1:k-1

        dec = [(sum(Linv[:,i].^2)/Linv[i,i]) for i = 1:n-dep]
        u = nodes[dep + 1]
        cf[dep + 1] = dec[u]
        Linv = delnode(Linv,u,n-dep)
    end
    println("\t CFCC achieved by augtree: ",  n/trace(Linv))
    logw(w, "\t Cfcc values ", cf);
    return n/trace(Linv)
end
