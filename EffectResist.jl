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
        cf = calculateNodeDists(distances, G.n)
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



datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1



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
#    @time exact(G,w)
#    @time exact(G,0,w)
    @time approx(G,w)
    @time approx(G,0,w)
end
    logw(w, "-------------------------------------------------------- ")
exit(0)
n = 8
m = 20
Line = Components3Graph(n,m)
println(Line)
#exact_distance = exact(Line,1,w) # same as er_mppinv(Line)
exact(Line,w)
exact(Line,0,w)
approx(Line,w)
approx(Line,0,w)
# LightGraphs structure is needed only to locate bridges
SLine =  sparseAdja2(Line)
bg = bridges(LightGraphs.Graph(SLine))
println(bg)
dist = zeros(Float64,n,n)
for e in bg
    x = e.src
    y = e.dst
    println(x,",",y)
    #pull!(nbr[x][y])
    SLine[x,y] = 0
    SLine[y,x] = 0
    dist[x,y] = 1
    dist[y,x] = 1
end
dropzeros!(SLine)
comps, nodes = allComp(SLine)
#println(comps)
println(nodes)
println("Running components:")
for c in comps
    if c.n != 1
        # ------- exact -------
        # println(c)
        C = sparsemat2Graph(c)
        ldist = exact(C,0,w)
        # ---------------------
        # ------- approx -------
        ### C = convert(SparseMatrixCSC{Float64,Int64},c)
        # println(sparsemat2Graph(c))
        # ldist = approx(c, w)
        # ---------------------
        #println("Cf values for components")
        #println(ldist)
        
    end
end
    
# dropzero! is needed to disconnect the components for real
#     println("after deletion of bridges")
#     println(full(spPath))
#     println(spPath)
#     comps, nodes = allComp(spPath)
#     println(nodes)
#     println("components")
#     for c in comps
#         if c.n != 1
#             local_er = approx(convert(SparseMatrixCSC{Float64,Int64},c),w)
#             #global_er   
#         end
#     end

    

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
