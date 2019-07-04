
using LightGraphs

mutable struct Component
    A:: SparseMatrixCSC{Float64}
    nc::Int64
    nodemap::Array{Int64,1} # size of nc Array{Int64,1}
    bdry::Array{Int64,1} # bdry nodes in Cluster
    link::Array{Int64,1}
    linkc::Int64
    distances:: Array{Float64}
    external::Array{Int64,1}
end

Component(A::SparseMatrixCSC{Float64},nodemap::Array{Int64,1}) = Component(A, A.n, nodemap,
                                                                           nothing, 0,
                                                                           #zeros(A.n)
                                                                           )

function addExtnodesInComponent(c:: Component, external:: Array{Int64,1})
    if c.external == nothing
        c.external = external
    end
end

# struct SimpleEdge{T<:Integer} #<: AbstractSimpleEdge{T}
#      src::T
#      dst::T
# end

struct Bridges
    #    edges:: Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}
    edges:: Array{Int64,1}
    m :: Int64 # number of edges in Bridges = size(edges)/2
    #    core1nodes :: Set{Int64} # I don't need to store core1nodes here
    #core2nodes :: Set{Int64}
    core2nodes :: Vector{Int64}
    n :: Int64 # number of core2nodes
    ext ::Array{Int64,1} # size of n , default 0, count of ext core1 nodes for each core2node 
    comp ::Array{Int64,1} # size of edges , corresponding component
end

Bridges(edges, core2nodes:: Set{Int64}, ext:: Array{Int64,1}) = Bridges(edges, length(edges)/2, core2nodes, length(core2nodes), ext, zeros(Int64,length(edges)))

Bridges(edges :: Array{Int64,1}, core2nodes:: Array{Int64,1}, ext:: Array{Int64,1}) = Bridges(edges, length(edges)/2, core2nodes, length(core2nodes), ext, zeros(Int64,length(edges)))


# Bridges(edges, nodes) =
#     Bridges(edges, size(edges,1),
#             nodes,
#             size(nodes,1),
#             zeros(Int64,size(nodes,1))
#             )

# Bridges(edges, nodes:: Set{Int64}) =
#     Bridges(edges, size(edges,1),
#             nodes,length(nodes),
#             zeros(Int64,length(nodes))
#             )


function printComponent(C:: Component)
    println("- A type:",eltype(C.A))
    println("- A.n=",C.nc)
    println("- list of nodes=",C.nodemap)
    println("- list of bdry=",C.bdry)
    println("- list of link=",C.link)
    println("- link.n=",C.linkc," (",100*C.linkc/C.nc,"%)")
    println("- list of external=",C.external)
    #println("- cf-distances=",C.distances)
end

function printEdges(edges:: Array{Int64,1}, n::Int64)
    print("- edges: ")
    for i in 1:2:n
        print("(",edges[i],",",edges[i+1],") ")
    end
    println("")
end


function printBridges(B:: Bridges)
    #println()
    printEdges(B.edges,length(B.edges))
    println("- m=",B.m)
    println("- length of edges =",length(B.edges))
    println("- n=",B.n)
    println("- list of core2nodes=",B.core2nodes)
    println("- list of ext (count)=",B.ext)
    println("- list of core3nodes=",B.edges)
    println("- comp of each node=",B.comp)
end

