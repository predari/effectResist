using LightGraphs
include("structures.jl")
"""
    bridges(g)
Compute the [bridges](https://en.m.wikipedia.org/wiki/Bridge_(graph_theory))
of a connected graph `g` and return an array containing all bridges, i.e edges
whose deletion increases the number of connected components of the graph.
# Examples
```jldoctest
julia> using LightGraphs
julia> bridges(StarGraph(5))
8-element Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}:
 Edge 1 => 2
 Edge 1 => 3
 Edge 1 => 4
 Edge 1 => 5
julia> bridges(PathGraph(5))
8-element Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}:
 Edge 4 => 5
 Edge 3 => 4
 Edge 2 => 3
 Edge 1 => 2
```
"""




function bridges end
#@traitfn
#function bridges(g::AG::(!IsDirected)) where {T, AG<:AbstractGraph{T}}
function bridges(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    bridges = Edge{T}[]   #keeps record of the bridge-edges
    
    # We iterate over all vertices, and if they have already been visited (pre != 0), we don't start a DFS from that vertex.
    # The purpose is to create a DFS forest.
    @inbounds for u in vertices(g)
        pre[u] != 0 && continue
        v = u #currently visiting vertex
        wi::T = zero(T) #index of children of v
        w::T = zero(T) #children of v
        cnt::T = one(T) # keeps record of the time
        first_time = true
        
        #start of DFS
        while !isempty(s) || first_time
            first_time = false
            if  wi < 1 #initialisation for vertex v
                pre[v] = cnt
                cnt += 1
                low[v] = pre[v]
                v_neighbors = outneighbors(g, v)
                wi = 1
            else
                wi, u, v = pop!(s) # the stack states, explained later
                v_neighbors = outneighbors(g, v)
                w = v_neighbors[wi]
                low[v] = min(low[v], low[w]) # condition check for (v, w) being a tree-edge
                if low[w] > pre[v]
                    edge = v < w ? Edge(v, w) : Edge(w, v)
                    push!(bridges, edge)
                end
                wi += 1
            end
            
            # here, we're iterating of all the childen of vertex v, if unvisited, we start a DFS from that child, else we update the low[v] as the edge is a back-edge.
            while wi <= length(v_neighbors)
                w = v_neighbors[wi]
                # If this is true , this indicates the vertex is still unvisited, then we push this on the stack.
                # Pushing onto the stack is analogous to visiting the vertex and starting DFS from that vertex.
                if pre[w] == 0
                    push!(s, (wi, u, v)) # the stack states are (index of child, currently visiting vertex, parent vertex of the child)
                    #updates the value for stimulating DFS from top of the stack
                    wi = 0 
                    u = v
                    v = w
                    break
                elseif w != u # (v, w) is a back-edge
                    low[v] = min(low[v], pre[w]) # condition for back-edges
                end
                wi += 1
            end
            wi < 1 && continue
        end
        
    end
    
    return bridges
end
# second function with different return type

function bridges(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    ### if nodes of bridges belong to terminal components,
    ### we do not need to store the edge. Just the nodes
    core1nodes = Set{T}()
    core1neighbor = Array{Int64,1}()
    #core3nodes = Set{T}()
    
    # We iterate over all vertices, and if they have already been visited (pre != 0), we don't start a DFS from that vertex.
    # The purpose is to create a DFS forest.
    @inbounds for u in vertices(g)
        pre[u] != 0 && continue
        v = u #currently visiting vertex
        wi::T = zero(T) #index of children of v
        w::T = zero(T) #children of v
        cnt::T = one(T) # keeps record of the time
        first_time = true
        
        #start of DFS
        while !isempty(s) || first_time
            first_time = false
            if  wi < 1 #initialisation for vertex v
                pre[v] = cnt
                cnt += 1
                low[v] = pre[v]
                v_neighbors = outneighbors(g, v)
                wi = 1

            else
                wi, u, v = pop!(s) # the stack states, explained later
                v_neighbors = outneighbors(g, v)
                w = v_neighbors[wi]
                low[v] = min(low[v], low[w]) # condition check for (v, w) being a tree-edge
                # println("v=$v, w=$w ==>", v_neighbors, outneighbors(g, w))
                # println(pre[v] ," < ", low[w], " ==> bridge ")
                # if low[w] > pre[v]
                #     push!(edges,w)
                #     push!(edges,v)
                # end
                if low[w] > pre[v]
                    if length(v_neighbors) == 1
                        push!(core1nodes, v)
                        push!(core1neighbor, w)
                    elseif length(outneighbors(g, w)) == 1
                        push!(core1nodes, w)
                        push!(core1neighbor, v)
                    else
                        push!(edges,w)
                        push!(edges,v)
                    end
                end
                wi += 1
            end
            
            while wi <= length(v_neighbors)
                w = v_neighbors[wi]
                if pre[w] == 0
                    push!(s, (wi, u, v))
                    wi = 0 
                    u = v
                    v = w
                    break
                elseif w != u # (v, w) is a back-edge
                    low[v] = min(low[v], pre[w]) # condition for back-edges
                end
                wi += 1
            end
            wi < 1 && continue
        end
        
    end
    ## core2pairs is replaced by sizes
    # core2pairs = zeros(Int64, length(core2nodes))
    # for (idx, u) in enumerate(core2nodes)
    #     c = count(x->x==u,core1neighbor);
    #     core2pairs[idx] = c;
    # end
    # println(core1neighbor)
    # for i in core1neighbor
    #     if length(outneighbors(g, i)) == 2
    #         println(i," ",length(outneighbors(g, i))," " )
    #     end
    # end
    t = time()
    ###### TODO: sort(unique()) or unique(sort()) ?
    ### TODO: uncomment
    score1neighbor = sort(core1neighbor)
    core2nodes = Array{Int64}()
    core2nodes = unique(score1neighbor)
     sizes = zeros(Int64, length(core2nodes))
    idx = 1
    sizes[1] = 1
    for i in 2:length(score1neighbor)
        if score1neighbor[i - 1] == score1neighbor[i]
            sizes[idx] += 1
        else
            idx += 1
            sizes[idx] += 1
        end
    end
    ### TODO: sizes have been found but!! their index does not correspond
    ### to the index of core2nodes. ## One solution is to reorder sizes
    ### or core2nodes so that they much. The other solution is to sort
    ### everything (maybe use OrderedSets). Also, I do really need to store
    ### core2nodes and core1neighbors while searching for bridges... keep one and try to prealloc
    println("COUNT TIMING IS ",time()- t, "(s)")

    # println("core1nodes", core1nodes)
    # println("core2nodes", core2nodes)
    # println("core2pairs", core2pairs)
    # println("core3nodes", core3nodes)
    if length(findin(core1nodes, core2nodes)) != 0
        println("WARNING!! Problem with core2 calculation")
        exit(1)
    end
    if length(findin(core1nodes,edges)) != 0
        println("WARNING!! Problem with core3 calculation")
        exit(1)
    end
    #println(core2nodes)
    B = Bridges(edges, core2nodes, sizes)
    ### remove the creation of core3nodes!!! and call following constructor
    #B = Bridges(edges, core2nodes, core2pairs)
    
    return B, core1nodes
end

function nbofbridges(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    bridges = Edge{T}[]   #keeps record of the bridge-edges
    nbofBridges :: Int64 = 0
    
    @inbounds for u in vertices(g)
        pre[u] != 0 && continue
        v = u #currently visiting vertex
        wi::T = zero(T) #index of children of v
        w::T = zero(T) #children of v
        cnt::T = one(T) # keeps record of the time
        first_time = true
        
        #start of DFS
        while !isempty(s) || first_time
            first_time = false
            if  wi < 1 #initialisation for vertex v
                pre[v] = cnt
                cnt += 1
                low[v] = pre[v]
                v_neighbors = outneighbors(g, v)
                wi = 1
            else
                wi, u, v = pop!(s) # the stack states, explained later
                v_neighbors = outneighbors(g, v)
                w = v_neighbors[wi]
                low[v] = min(low[v], low[w]) # condition check for (v, w) being a tree-edge
                if low[w] > pre[v]
                    edge = v < w ? Edge(v, w) : Edge(w, v)
                    push!(bridges, edge)
                    nbofBridges += 1
                end
                wi += 1
            end
            
            # here, we're iterating of all the childen of vertex v, if unvisited, we start a DFS from that child, else we update the low[v] as the edge is a back-edge.
            while wi <= length(v_neighbors)
                w = v_neighbors[wi]
                # If this is true , this indicates the vertex is still unvisited, then we push this on the stack.
                # Pushing onto the stack is analogous to visiting the vertex and starting DFS from that vertex.
                if pre[w] == 0
                    push!(s, (wi, u, v)) # the stack states are (index of child, currently visiting vertex, parent vertex of the child)
                    #updates the value for stimulating DFS from top of the stack
                    wi = 0 
                    u = v
                    v = w
                    break
                elseif w != u # (v, w) is a back-edge
                    low[v] = min(low[v], pre[w]) # condition for back-edges
                end
                wi += 1
            end
            wi < 1 && continue
        end
        
    end
    
    return nbofBridges
end



function bridges2(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    ### if nodes of bridges belong to terminal components,
    ### we do not need to store the edge. Just the nodes
    ### core1nodes = Set{T}()
    ### core1neighbor = Array{Int64,1}()
    ### core3nodes = Set{T}()
    m :: Int64 = 0
    tms = zeros(T, nv(g))
    degree1neighbor = Array{Int64,1}()


    @inbounds for u in vertices(g)
        pre[u] != 0 && continue
        v = u #currently visiting vertex
        wi::T = zero(T) #index of children of v
        w::T = zero(T) #children of v
        cnt::T = one(T) # keeps record of the time
        first_time = true
        
        #start of DFS
        while !isempty(s) || first_time
            first_time = false
            if  wi < 1 #initialisation for vertex v
                pre[v] = cnt
                cnt += 1
                low[v] = pre[v]
                v_neighbors = outneighbors(g, v)
                wi = 1
            else
                wi, u, v = pop!(s) # the stack states, explained later
                v_neighbors = outneighbors(g, v)
                w = v_neighbors[wi]
                low[v] = min(low[v], low[w]) # condition check for (v, w) being a tree-edge
                #println("v=$v, w=$w ==>", v_neighbors, outneighbors(g, w))
                #println(pre[v] ," < ", low[w], " ==> bridge ")
                if low[w] > pre[v]
                    push!(edges,w)
                    tms[w] += 1
                    push!(edges,v)
                    tms[v] += 1
                    if length(v_neighbors) == 1
                        push!(degree1neighbor, w)
                    elseif length(outneighbors(g, w)) == 1
                        push!(degree1neighbor, v)
                    end

                end
                wi += 1
            end

            
            while wi <= length(v_neighbors)
                w = v_neighbors[wi]
                if pre[w] == 0
                    push!(s, (wi, u, v))
                    wi = 0 
                    u = v
                    v = w
                    break
                elseif w != u # (v, w) is a back-edge
                    low[v] = min(low[v], pre[w]) # condition for back-edges
                end
                wi += 1
            end
            wi < 1 && continue
        end
        
    end
    t = time()
    #core1nodes = Array{Int64,1}()
    core1nodes = Set{T}()
    links = Array{Int64,1}()
    linksneighbor = Array{Int64,1}()
    uedges = unique(edges)
    m = length(edges)
    core2nodes = Array{Int64}()
    core2nodes = unique(degree1neighbor)
    cntrcore2nodes = tms[core2nodes]
    remove = Array{Int64,1}()

    for i in uedges
        tms[i] = length(outneighbors(g, i)) - tms[i] + 1  
    end
    for (idx, u) in enumerate(core2nodes)
        if tms[u] == 1
            push!(remove,idx)
        end
    end
    deleteat!(core2nodes,remove)
    nbc2 :: Int64 = length(core2nodes) 
    deleteat!(cntrcore2nodes,remove)
    
    for i in 1:2:m
        if tms[edges[i]] == 1 && tms[edges[i + 1]] == 1
            push!(core1nodes,edges[i])
            push!(core1nodes,edges[i+1])
        elseif tms[edges[i]] > 1 && length(outneighbors(g, edges[i+1])) > 1
            # category of rdg(i,i+1) = (C3,C2) falls also here
            push!(links, edges[i])
            push!(linksneighbor, edges[i+1])
        elseif tms[edges[i + 1]] > 1  && length(outneighbors(g, edges[i])) > 1
            push!(links, edges[i + 1])
            push!(linksneighbor, edges[i])
        end
    end
    ### setdiff(u,v): here u must be the smaller set
    ### we need space for nodes that are links with
    ### paths and not nodes of type core2nodes
    onlylinks = length(setdiff(links, core2nodes))
    sizes = Array{Array{Int, 1}}(nbc2 + onlylinks)
    for i in indices(sizes,1) sizes[i] = []; end
    for i in 1 : nbc2
        push!(sizes[i], cntrcore2nodes[i])
    end
    println(sizes)
    

    duplicate = zeros(Int64,length(linksneighbor))    
    for (idx, u) in enumerate(unique(linksneighbor))
        for v in linksneighbor
            if u == v
                if duplicate[idx] > 0
                    tms[u] = 2 ## 2 is a random number. Could be any more than 1
                    break;
                end
                duplicate[idx] += 1
            end
        end
    end

bridges =  Array{Int64,1}()
i :: Int64 = 1
for (idx,l) in enumerate(links)
    if tms[linksneighbor[idx]] > 1
        push!(bridges,l)
        push!(bridges,linksneighbor[idx])
    else
        p,d = bfs_edge_subtree2(g, l, linksneighbor[idx], tms)        
        if length(d) == 1
            ## TODO: assert length(path) == 2
            push!(bridges,p[1])
            push!(bridges,p[2])
        else 
            cntr = zeros(Int64, maximum(d))
            cntr = count2(d,length(d))
            if isempty(findin(core2nodes, l))
                push!(sizes[nbc2 + i], cntr)
                i += 1
            else 
                idx2 = getindex(findin(core2nodes, l))
                for j in 2:length(cntr)
                    push!(sizes[idx2],cntr[j])
                end
            end
        end
    end
end

println("core2nodes:", core2nodes)
println("sizes:", sizes)
for i in 1:length(sizes)
    println(sizes[i])
end
println("final bridges = ", bridges)
B = Bridges(bridges, core2nodes, sizes)
println(length(core1nodes))
#return B, core1nodes
shell1nodes = Array{Int64,1}
shell1nodes = k_shell(g, 1)
println(length(shell1nodes))
return B, shell1nodes
end

### TODO start BFS from linksneighbor
### TODO stop if node has tms > 1 or if length(outneighbors(g, v) == 1

function count2( array ::  Array{Int64,1}, len :: Int64)
    cntr = zeros(Int64, maximum(array))
    #cntr =  Array{Int64,1}()
    c :: Int64 = 1 
    for (idx,u) in enumerate(array)
        if idx == len
            cntr[u] = c
        elseif idx > 1
            if u == array[idx - 1]
                c += 1
            else
                cntr[array[idx - 1]] = c
                c = 1
            end
        end
    end
    return cntr
end

function bfs_edge_subtree2(g::AbstractGraph{T}, source :: Int64, next:: Int64, tms ::  Array{Int64,1} ) where T
    
    n = nv(g)
    visited = falses(n)
    distance = Vector{T}()
    path = Vector{T}()
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    visited[source] = true
    visited[next] = true
    push!(cur_level, next)
    cnt :: Int64 = 1
    push!(distance,cnt)
    push!(path,next)
    cnt += 1
   
    while !isempty(cur_level)
        @inbounds for v in cur_level
# check what @simd is for
            @inbounds @simd for i in  outneighbors(g, v)
                if !visited[i] && tms[i] == 1
                    push!(next_level, i)
                    push!(distance,cnt)
                    push!(path,i)
                    visited[i] = true
                elseif !visited[i] && tms[i] != 1
                    ## this code is for paths that are connected to
                    ## a component of 1 node, since this node is also
                    ## connected to other components. Example:
                    ##        C1 -- 2 -- C3
                    ##              |
                    ##              4
                    ##              |
                    ##             C6
                    ## the solution is to consider 2--4 a bridge and
                    ## 4 is added to component C6 without considering
                    ## edge 4--6 as a bridge. This is done to make things
                    ## simpler as I don't think that they will impact
                    ## a lot the algorithm.
                    println("WARNING: code may has some logical bugs!")
                    push!(distance,cnt)
                    last = path[end]
                    empty!(path)
                    push!(path,last)
                    push!(path,i)
                    return path , 1
                end
            end
        end
        empty!(cur_level)
        cnt += 1
        cur_level, next_level = next_level, cur_level
        ## do I need sorting?
        #sort!(cur_level)
    end
    return path,distance
end


function bfs_edge_subtree(g::AbstractGraph{T}, source :: Int64, next:: Int64) where T
    n = nv(g)
    visited = falses(n)
    distance = Vector{T}()
    path = Vector{T}()
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    visited[source] = true
    visited[next] = true
    push!(cur_level, next)
    cnt :: Int64 = 1
    push!(distance,cnt)
    push!(path,next)
    cnt += 1
   
    while !isempty(cur_level)
        @inbounds for v in cur_level
            @inbounds @simd for i in  outneighbors(g, v)
                if !visited[i]
                    push!(next_level, i)
                    push!(distance,cnt)
                    push!(path,i)
                    visited[i] = true
                end
            end
        end
        empty!(cur_level)
        cnt += 1
        cur_level, next_level = next_level, cur_level
        ## do I need sort
        sort!(cur_level)
    end
    println("path $source:",path)
    println("distance $source:",distance)
    return path,distance    
end

bfs_parents2(g::AbstractGraph, s::Integer; dir = :out) = 
    (dir == :out) ? _bfs_parents(g, s, outneighbors) : _bfs_parents(g, s, inneighbors)

function _bfs_parents(g::AbstractGraph{T}, source, neighborfn::Function) where T
    n = nv(g)
    visited = falses(n)
    parents = zeros(T, nv(g))
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    @inbounds for s in source
        visited[s] = true
        push!(cur_level, s)
        parents[s] = s
    end
    while !isempty(cur_level)
        @inbounds for v in cur_level
            @inbounds @simd for i in  neighborfn(g, v)
                if !visited[i]
                    push!(next_level, i)
                    parents[i] = v
                    visited[i] = true
                end
            end
        end
        empty!(cur_level)
        cur_level, next_level = next_level, cur_level
        sort!(cur_level)
    end
    return parents
end
