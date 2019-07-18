include("graph.jl")
include("lapl.jl")
include("logw.jl")

include("fastapprox.jl")
include("methods.jl")
include("sampling.jl")
using Laplacians


datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1
logw(w, "-------------------------------------------------------- ")
# #n = 8
# #m = 20
# #Line = Components3Graph(n,m)

# #Line = TestGraph(20, 54)
# #Line = TestGraph(18, 48)
# Line = TestGraph(21, 54)
# println(Line)
# A, L = sparseAdja(Line)
# #@time approx(A,L,w)
# #Line = Components3Graph(8, 22)
# #Line = TestGraph(21, 54)
# #Line = TestGraph(18, 48)
# #Line = TestGraph(20, 54)
# #println(Line)
# #A, L = sparseAdja(Line)

# #@time max = exact(Line,w)

Line = TestGraph(21, 54)
println(Line)
A, L = sparseAdja(Line)
@time cfcAccelerate(A, w, 25)

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

    A, L = sparseAdja(G)
#    @time  max = approx(A,L,w)
#    @time max = exact(G,w)
#    @time cfcAccelerate(A, w, 25)    
#    @time approxcore2(A, L, w)
#    A, L = sparseAdja(G)
    @time cfcAccelerate(A, w, 25)    
end

# - list of core2nodes=[1, 211, 289, 290, 999, 1000, 1135, 2134, 2147, 2792]
# - list of ext (count)=[280, 455, 756, 706, 170, 92, 57, 31, 147, 96]
# - list of core3nodes=[2134, 2075, 2075, 289, 1135, 1020, 1020, 1000, 2792, 2384, 2384, 211]
# - comp of each node=[2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]

# - list of core2nodes=[1, 289, 2134, 2147, 1000, 999, 211, 290]
# - list of ext (count)=Array{Int64,1}[[280], [757], [32], [147], [93, 1, 56], [170], [456, 1, 95], [706]]
# - list of core3nodes=[2134, 2075, 289, 2075]
