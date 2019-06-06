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
    println("Distances:", distances)
    if alldistances == 1
        return distances
    else
        if size(distances,1) == G.n-1
            cf = calculateNodeDists(distances, G.n)
        else
            cf = distances
        end
        #println("cf size:", size(cf,1))
        cf = calculateCF(cf, G.n)
        # cf = zeros(G.n)
        # for i in 1:G.n
        #     cf[i] = G.n/distances[i]
        # end
        #println(cf)
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
        println("Distances:", cf)
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
#     @time exact(G,w)
#     @time exact(G,0,w)
#     @time approx(G,w)
#     @time approx(G,0,w)
# end
# logw(w, "-------------------------------------------------------- ")

# L = LineGraph(7)
# println(L)
# distances = exact_comp(L,w)
n = 8
m = 20
Line = Components3Graph(n,m)
println(Line)
#exact_distance = exact(Line,1,w) # same as er_mppinv(Line)
#exact(Line,w)
exact(Line,0,w)
#approx(Line,w)
approx(Line,0,w)
#distances = exact_comp(Line, w)
#println(distances)
nbr = Array{Array{Int, 1}}(4)
for i in indices(nbr,1) nbr[i] = [] end
push!(nbr[1], 2)
push!(nbr[1], 4)
push!(nbr[2], 1)
push!(nbr[2], 4)
push!(nbr[2], 3)
push!(nbr[3], 2)
push!(nbr[3], 4)
push!(nbr[4], 1)
push!(nbr[4], 3)
push!(nbr[4], 2)
Line = Graph(4,10,nbr)
println(Line)
#exact(Line,w)
exact(Line,0,w)
#approx(Line,w)
approx(Line,0,w)
Line = TriangleGraph(3,6)
println(Line)
#exact(Line,w)
exact(Line,0,w)
#approx(Line,w)
approx(Line,0,w)

# Line = ComponentExtnodes3Graph(7, 16)
# println(Line)
# exact(Line,w)
# exact(Line,0,w)
# approx(Line,w)
# approx(Line,0,w)
# nbr = Array{Array{Int, 1}}(4)
# for i in indices(nbr,1) nbr[i] = [] end
# push!(nbr[1], 2)
# push!(nbr[1], 3)
# push!(nbr[2], 1)
# push!(nbr[2], 4)
# push!(nbr[2], 3)
# push!(nbr[3], 1)
# push!(nbr[3], 2)
# push!(nbr[3], 4)
# push!(nbr[4], 3)
# push!(nbr[4], 2)
# rLine = Graph(4,10,nbr)
# println(rLine)
# exact(rLine,w)
# exact(rLine,0,w)
# approx(rLine,w)
# approx(rLine,0,w)
# Line = Component2sExtnodes3Graph(10, 24)
# println(Line)
# exact(Line,w)
# exact(Line,0,w)
# approx(Line,w)
# approx(Line,0,w)


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
