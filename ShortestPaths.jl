include("structures.jl")
using Laplacians

mutable struct intHeap{Tkey,Tind}
  keys::Vector{Tkey}
  heap::Vector{Tind}
  index::Vector{Tind}
  nitems::Tind
end #intHeap

intHeap(n::Int64) = intHeap(Inf*ones(Float64,n),-ones(Int64,n),zeros(Int64,n),0)
intHeap(n::Int32) = intHeap(Inf*ones(Float32,n),-ones(Int32,n),zeros(Int32,n),0)

function intHeapAdd!(nh::intHeap, node::Tind, key::Tkey) where {Tkey,Tind}
  if nh.index[node] > 0 # if already in the heap
    if key < nh.keys[node]
      intHeapSet!(nh, node, key)
    end

  else # if it really is new

    nhp = nh.nitems+1

    nh.keys[node] = key
    nh.heap[nhp] = node
    nh.index[node] = nhp
    nh.nitems = nhp

    intHeapUp!(nh, node)

  end
end # intHeapAdd!

function intHeapAddGost!(nh::intHeap, node::Tind) where {Tind}
  if !(nh.index[node] > 0)
    nhp = nh.nitems+1

    nh.keys[node] = nh.keys[node-1]
    nh.heap[nhp] = node
    nh.index[node] = nhp
    nh.nitems = nhp
    intHeapUp!(nh, node)

  end
end # intHeapAdd!


function intHeapDown!(nh::intHeap, node::Tind) where Tind
  pos = nh.index[node]
  key = nh.keys[node]
  leftPos = pos*2
  moved = true
  @inbounds while (leftPos <= nh.nitems) && moved
    moved = false
    rightPos = pos*2+1

    if rightPos > nh.nitems
      childPos = leftPos
      childNode = nh.heap[childPos]
      childKey = nh.keys[childNode]
    else
      leftNode = nh.heap[leftPos]
      leftKey = nh.keys[leftNode]
      rightNode = nh.heap[rightPos]
      rightKey = nh.keys[rightNode]

      if leftKey < rightKey
        childPos = leftPos
        childNode = leftNode
        childKey = leftKey
      else
        childPos = rightPos
        childNode = rightNode
        childKey = rightKey
      end
    end

    if childKey < key
      nh.heap[childPos] = node
      nh.heap[pos] = childNode
      nh.index[node] = childPos
      nh.index[childNode] = pos

      pos = childPos
      leftPos = pos*2
      moved = true
    end

  end #while
end # intHeapDown!

function intHeapPop!(nh::intHeap)
  minNode = nh.heap[1]

  nh.index[minNode] = 0

  @inbounds if (nh.nitems > 1)
    node = nh.heap[nh.nitems]
    nh.heap[1] = node
    nh.index[node] = 1
    intHeapDown!(nh, node)
  end
  nh.nitems = nh.nitems - 1

  return minNode
end # intHeapPop!

function intHeapUp!(nh::intHeap, node::Tind) where Tind
  pos = nh.index[node]
  moved = true

  @inbounds while (pos > 1) && moved
    key = nh.keys[node]

    parentPos = div(pos,2)
    parentNode = nh.heap[parentPos]
    parentKey = nh.keys[parentNode]

    moved = false

    if (parentKey > key)
      nh.heap[parentPos] = node
      nh.heap[pos] = parentNode
      nh.index[node] = parentPos
      nh.index[parentNode] = pos
      pos = parentPos
      moved = true
    end
  end

end # intHeapUp!

function intHeapSort(x::Array{Float64})
  n = length(x)
  nh = intHeap(n)

  @inbounds for i in 1:n
    intHeapAdd!(nh, i, x[i])
  end

  out = zeros(Float64,n)
  @inbounds for i in 1:n
    out[i] = nh.keys[intHeapPop!(nh)]
  end

  return out

end # intHeapSort


function intHeapSort(nh::intHeap)
  n = length(nh.keys)

  out = zeros(Float64,n)
  for i in 1:n
    out[i] = nh.keys[intHeapPop!(nh)]
  end

  return out

end # intHeapSort

function intHeapSet!(nh::intHeap, node::Tind, key::Tkey) where {Tkey,Tind}
  oldKey = nh.keys[node]
  nh.keys[node] = key

  if (key < oldKey)
    intHeapUp!(nh,node)
  else
    intHeapDown!(nh,node)
  end
end # intHeapSet!

"""Computes the lenghts of shortest paths from `start`.
Returns both a vector of the lenghts, and the parent array
in the shortest path tree.
This algorithm treats edge weights as reciprocals of distances.
DOC BETTER
"""
function shortestPaths(mat::SparseMatrixCSC{Tv,Ti}, start::Ti, comp:: Array{Int64,1}, nc :: Int64) where {Tv,Ti}
    distances = zeros(Float64,nc)
    n = mat.n
    visited = zeros(Bool,n)
    visitedcomp = zeros(Bool,nc)

    nh = intHeap(n)
    dists = nh.keys
    #println("Dists",dists)
    pArray = zeros(Ti,n)

    intHeapAdd!(nh,start,0.0)
    pArray[start] = start

    while nh.nitems > 0
        v::Ti = intHeapPop!(nh)
        x = comp[v]
        println("vertex $v ... (comp $x)!")
        visited[v] = true
        visitedcomp[x] = true
        dv = dists[v]
#        println("visitedcomp = ", visitedcomp)
        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            nbr = mat.rowval[ind]
            y = comp[nbr]
            println("neighbor = $nbr (comp $y)")
            if !visited[nbr]
                #newdist = dv + 1/mat.nzval[ind]
                newdist = dv + mat.nzval[ind]
                if newdist < dists[nbr]
                    intHeapAdd!(nh,nbr,newdist)
                    if !visitedcomp[y]
                        distances[y] = newdist
                    end
                    visitedcomp[y] = true
                    #println("visitedcomp = ", visitedcomp)
                    #println(visited)
                    #println(dists)
                    pArray[nbr] = v
                end # if
            end # if
        end # for
        
    end # while
    println(dists)
    return distances, pArray
    #return copy(dists), pArray

end # shortestPaths


mutable struct cVertex
    id::Int64
    x::Array{Int64}
    rx::Array{Int64}
    multiplier::Array{Int64} # size of length(x)
    cmplist:: Array{Int64}
end 

cVertex(i::Int64,x::Array{Int64},rx::Array{Int64},n::Int64) = cVertex(i,x,rx,zeros(Int64,length(x)),zeros(Int64,n))
cVertex(i::Int64,n::Int64) = cVertex(i,[],[],zeros(Int64,1),zeros(Int64,n))

function printcVertex(cV :: cVertex)
    println("- cV.i=",cV.id)
    print("- list of x: ")
    println(cV.x)
    print("- list of rx: ")
    println(cV.rx)
    println("- multipliers: ")
    println(cV.multiplier)
    println("- cmplist: ")
    println(cV.cmplist)
end

# nodes -> newedges
function shortestPaths(mat::SparseMatrixCSC{Tv,Ti}, start::Ti, comp:: Array{Int64,1}, nodes:: Array{Int64,1}, nc :: Int64) where {Tv,Ti}
    # distances and path of components
    distances = zeros(Float64,nc)
    path = zeros(Int64,nc)
    path[comp[start]] = -1
    n = mat.n
    visited = zeros(Bool,n)
    visitedcomp = zeros(Bool,nc)
    nh = intHeap(n)
    dists = nh.keys

    onetime = true
    pArray = zeros(Ti,n)

    intHeapAdd!(nh,start,0.0)
    pArray[start] = start
    
    while nh.nitems > 0
        v::Ti = intHeapPop!(nh)
        x = comp[v]
        println("vertex $v ... (comp $x)")
        visited[v] = true
        if path[x] == 0
            path[x] = nodes[v]
        end
        visitedcomp[x] = true
        dv = dists[v]
 #       println("visitedcomp = ", visitedcomp)
        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            nbr = mat.rowval[ind]
            y = comp[nbr]
           println("neighbor = $nbr (comp $y)")
            if !visited[nbr]
                #newdist = dv + 1/mat.nzval[ind]
                newdist = dv + mat.nzval[ind]
                if newdist < dists[nbr]
                    intHeapAdd!(nh,nbr,newdist)
                    if !visitedcomp[y]
                        distances[y] = newdist
                    end
                    visitedcomp[y] = true
                    #println("visitedcomp = ", visitedcomp)
                    #println(visited)
                    #println(dists)
                    pArray[nbr] = v
                end # if
            end # if
        end # for
        
    end # while
    println(path)
    return distances, path

end # shortestPaths
