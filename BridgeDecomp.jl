include("graph.jl")
include("lapl.jl")
include("logw.jl")

include("fastapprox.jl")
include("methods.jl")
include("sampling.jl")
using Laplacians

function callCore2Approx(G :: Graph, w :: IOStream, solution = nothing)
    A = sparseAdja2(G)
    wrapproxcore2(A, w, solution)
end

function callSamplingApprox(G :: Graph, w :: IOStream, pv:: Int64, solution = nothing)
    A = sparseAdja2(G)
    wrsamplingApprox(A, w, pv, solution)
end

function callApprox(G :: Graph, w :: IOStream, solution = nothing)
    A, L = sparseAdja(G)
    wrapprox(A,L,w, solution)
end

function callFastApprox(G :: Graph, w :: IOStream, solution = nothing)
    A = sparseAdja2(G)
    wrcfcAccelerate(A, w, solution)
end

datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")

e = length(ARGS) >= 2 ? parse(Float64,ARGS[2]) : 0.1
logw(w, "-------------------------------------------------------- ")

#Line = Components3Graph(8, 22)
#Line = TestGraph(18, 48)
Line = TestGraph(21, 54)
println(Line)
#@time max = exact(Line,w)
A, L = sparseAdja(Line)
wrapprox(A,L,w)
A, L = sparseAdja(Line)
wrcfcAccelerate(A, w)
A, L = sparseAdja(Line)
wrapproxcore2(A, w)

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
#    callApprox(G, w)
#    callSamplingApprox(G, w, 10)
    @time    callCore2Approx(G, w)
    @time    callFastApprox(G, w)
   exit() 
end

# - list of core2nodes=[1, 211, 289, 290, 999, 1000, 1135, 2134, 2147, 2792]
# - list of ext (count)=[280, 455, 756, 706, 170, 92, 57, 31, 147, 96]
# - list of core3nodes=[2134, 2075, 2075, 289, 1135, 1020, 1020, 1000, 2792, 2384, 2384, 211]
# - comp of each node=[2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]

# - list of core2nodes=[1, 289, 2134, 2147, 1000, 999, 211, 290]
# - list of ext (count)=Array{Int64,1}[[280], [757], [32], [147], [93, 1, 56], [170], [456, 1, 95], [706]]
# - list of core3nodes=[2134, 2075, 289, 2075]
