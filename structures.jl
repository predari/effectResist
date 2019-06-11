
using LightGraphs

mutable struct Component
    A:: SparseMatrixCSC{Float64}
    nc::Int64
    nodemap::Array{Int64,1} # size of nc Array{Int64,1}
    bdry::Array{Int64,1} # bdry nodes in Cluster
    bdryc::Int64
    distances:: Array{Float64}
    external::Array{Int64,1}
end

# struct SimpleEdge{T<:Integer} #<: AbstractSimpleEdge{T}
#      src::T
#      dst::T
# end

struct Bridges
#    edges :: Array{SimpleEdge{Int64},1} #bridges
    edges:: Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}
    m :: Int64 # number of edges in Bridges
    nodes :: Set{Int64} # nodes related to Bridges instead of Array{Int64,1}
    n :: Int64 # number of nodes
    comp ::Array{Int64,1} # size of n , corresponding component
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


                                                                         
Bridges(edges, nodes) =
    Bridges(edges, size(edges,1),
            nodes,
            size(nodes,1),
            zeros(Int64,size(nodes,1))
            )

Bridges(edges, nodes:: Set{Int64}) =
    Bridges(edges, size(edges,1),
            nodes,length(nodes),
            zeros(Int64,length(nodes))
            )

function printComponent(C:: Component)
    println("- A type:",eltype(C.A))
    println(C.A)
    println("- A.n=",C.nc)
    println("- list of nodes=",C.nodemap)
    println("- list of bdry=",C.bdry)
    println("- bdry.n=",C.bdryc," (",100*C.bdryc/C.nc,"%)")
    println("- list of external=",C.external)
end

function printEdges(edges:: Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1})
    m = size(edges,1)
    
    print("- edges: ")
    for i = 1:m
        print("(",edges[i].src,",",edges[i].dst,") ")
    end
    println("")
end


function printBridges(B:: Bridges)
    #println()
    printEdges(B.edges)
    println("- m=",B.m)
    println("- list of nodes=",B.nodes)
    println("- n=",B.n)
    println("- comp of each node=",B.comp)
end

