include("graph.jl")
include("lapl.jl")
include("logw.jl")
include("bridge.jl")
include("components.jl")
include("Exact.jl")
include("appxInvTrace.jl")
include("Linvdiag.jl")
include("Lprimitives.jl")
include("compDistances.jl")
using Laplacians
using LightGraphs

function isTree(A::SparseMatrixCSC)
    isConnected(A) && (nnz(A) == 2*(A.n-1))
end


function approx(G, alldistances, w :: IOStream)
    logw(w,"****** Running approx ******")
    distances = erJLT(G, 1)
    if alldistances == 1
        return distances
    else
        if size(distances,1) == G.n-1
            cf = calculateNodeDists(distances, G.n)
        else
            cf = distances
        end
        println("cf size:", size(cf,1))
        cf = calculateCF(cf, G.n)
        # cf = zeros(G.n)
        # for i in 1:G.n
        #     cf[i] = G.n/distances[i]
        # end
        logw(w,"\t node with argmax{c(", indmax(cf), ")} = ",
             maximum(cf))
        return cf
    end
end

function approx(G, w :: IOStream)
    logw(w,"****** Running (chinese) approx ******")
    #start_time = time()
    u, maxcf = erJLT(G)
    #elapsed_time = time()-start_time
    logw(w,"\t node with argmax{c(", u, ")} = ", maxcf)
    #logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")
end




function exact(G, alldistances, w :: IOStream)
    logw(w,"****** Running (my) exact ******")
    # time start
    # start_time = time()
    distances = erINV(G,1)
    if alldistances == 1
        # elapsed_time = time()-start_time
        # logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")
        return distances
    else
        cf = calculateNodeDists(distances, G.n)
        #println(cf)
        cf = calculateCF(cf, G.n)
        logw(w,"\t node with argmax{c(", indmax(cf), ")} = ",
             maximum(cf))
        #elapsed_time = time()-start_time
        #logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")
        return cf
    end
end

function exact(G, w :: IOStream)
    logw(w,"****** Running (chinese) exact ******")
    # time start
    #start_time = time()
    u, maxcf = erINV(G)
    #elapsed_time = time()-start_time
    logw(w,"\t node with argmax{c(", u, ")} = ", maxcf)
    #logw(w,"\t totaltime: ",time() - elapsed_time, " (s)")

end


function exact_comp(G, w :: IOStream)
    logw(w,"****** Running exact components ******")
    A =  sparseAdja2(G)
    n = G.n
    bg = bridges(LightGraphs.Graph(A))
    println(bg)
    distances = zeros(n,n)
    # distances = Array{Array{Float64, 1}}(n-1)
    # for i in indices(distances,1) distances[i] = [] end
    # for i in 1:n-1
    #     for j in i+1:n
    #         push!(distances[i], 0.0)
    #     end
    # end
    for e in bg
        x = e.src
        y = e.dst
        #pull!(nbr[x][y])
        A[x,y] = 0
        A[y,x] = 0
    end
    dropzeros!(A)
    comps, nodes = allComp(A)
    println(comps)
    println(nodes)
    ncomps = 1
    #println("Running components:")
    for c in comps
        println(c,ncomps)
        if c.n != 1
            C = sparsemat2Graph(c)
            idx = nodes[ncomps]
            println(idx)
            ldistance = exact(C,1,w)
            ldistance = full(ldistance)
            s = size(ldistance,1)
            println("component size:",s, ",", size(ldistance,2))
            println(ldistance)
            for i in 1:s
                distances[idx[i],idx[i+1:end]] = ldistance[i]
                println(idx[i],",",idx[i+1:end])
                println(ldistance[i])
            end
        end
        ncomps = ncomps + 1
    end
    for i in 1:n
        for j in 1:n
            if distances[i,j] != 0.0
                distances[j,i] = distances[i,j]
            end
        end
    end
    
    println(ncomps)
    idx = []
    for e in bg
        x = e.src
        y = e.dst
        println(x,",",y)
        #pull!(nbr[x][y])
        if x < y 
            distances[x,y] = 1
            distances[y,x] = 1
            for i in 1:ncomps
                j = find(e -> e == x ,nodes[i])
                if isempty(j) == false
                    
                    c = comps[i]
                    idx = nodes[i]
                    #println("(",c,",",idx,")")
                    break; # (i,j) (component, idx in nodes)
                end
            end
            #println("(",i,",",j,")")
            idxr = filter!(e->e!=x,idx)
            println(idxr)
            for l in idxr
                distances[l,y] = distances[l,x] + 1
            end
            for i in 1:n
                for j in 1:n
                    if distances[i,j] != 0.0
                        distances[j,i] = distances[i,j]
                    end
                end
            end
        end
    end
    distances
return distances
end


datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1



# for rFile in filter( x->!startswith(x, "."), readdir(string(datadir)))
#     logw(w, "---------------------",rFile,"-------------------------- ")
#     logw(w, "Reading graph from edges list file ", rFile)
#     G = read_file(string(datadir,rFile))
#     logw(w, "\t G.n: ", G.n, "\t G.m: ", G.m)
#     A, L = sparseAdja(G)
#     if !Laplacians.isConnected(A)
#         logw(w," WARNING: Graph is not connected. Program will exit!");
#         exit()
#     end
# #    @time exact(G,w)
# #    @time exact(G,0,w)
# #    @time approx(G,w)
#     @time approx(G,0,w)
# end
    logw(w, "-------------------------------------------------------- ")

# L = LineGraph(7)
# println(L)
# distances = exact_comp(L,w)
n = 8
m = 20
Line = Components3Graph(n,m)
println(Line)
#exact_distance = exact(Line,1,w) # same as er_mppinv(Line)
#exact(Line,w)
#exact(Line,0,w)
#approx(Line,w)
#approx(Line,0,w)
distances = exact_comp(Line, w)
println(distances)


    #T1 = akpwU(A);
    #if !isTree(T1)
    #    logw(w,"WARNING: Graph is a tree! ");
    #end
    #comp = randperm(MersenneTwister(1234), A.n)
    #compGraphU(A, comp)
    logw(w,"****** Running starDecomp ******")
    #T = starDecomp(A)
    #print(T)
    # if !isTree(T)
    #     logw(w,"WARNING: Graph is a tree! ");
    # end
    # logw(w, "\t T.n: ", T.n, "\t T.m: ", nnz(T))

close(w)
