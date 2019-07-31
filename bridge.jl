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

###### ORIGINAL BRIDGES FUNCTION FROM LIGHTGRAPHS #######
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

### old bridges function that used to calculate only core 1 connections
### to bndr nodes.
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
    ### core2nodes and core1neighbors while searching for bridges...
    ### keep one and try to prealloc
    println("COUNT TIMING IS ",time()- t, "(s)")

    if length(findin(core1nodes, core2nodes)) != 0
        println("WARNING!! Problem with core2 calculation")
        exit(1)
    end
    if length(findin(core1nodes,edges)) != 0
        println("WARNING!! Problem with core3 calculation")
        exit(1)
    end
    B = Bridges(edges, core2nodes, sizes)
    return B, core1nodes
end

### main bridges function that is used to calculate bridges
### the bndry nodes and all the related to bndry nodes
### simple paths of core1 type. 
function bridges2(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    m :: Int64 = 0
    tms = zeros(T, nv(g))
    degree1neighbor = Array{Int64,1}()

    st = time()
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
    
    println("### finding all bridges time:",time()- st, "(s)")
    t = time()
    core1nodes = Set{T}()
    links = Array{Int64,1}()
    linksneighbor = Array{Int64,1}()
    uedges = unique(edges)
    m = length(edges)
    core2nodes = Array{Int64}()
    core2nodes = unique(degree1neighbor)
    cntrcore2nodes = tms[core2nodes]
    remove = Array{Int64,1}()
    #println("core2nodes", core2nodes )
    #println("cntrcore2nodes", cntrcore2nodes )
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
    #println("After deleting some nodes")
    #println("core2nodes", core2nodes )
    #println("cntrcore2nodes", cntrcore2nodes )
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
    #println("links",links)
    #println("links length = ", length(links))
    #println("linksneighbor", linksneighbor)
    
    ### setdiff(u,v): here u must be the smaller set
    ### we need space for nodes that are links with
    ### paths and not nodes of type core2nodes
    onlylinks = length(setdiff(links, core2nodes))
    #println("onlylinks $onlylinks")
    #println(nbc2)
    sizes = Array{Array{Int, 1}}(nbc2 + onlylinks)
    for i in indices(sizes,1) sizes[i] = []; end
    for i in 1 : nbc2
        push!(sizes[i], cntrcore2nodes[i])
    end
    #println(sizes)
    #println(size(sizes,1))

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
    #println(length(links))
    #println(setdiff(links, core2nodes))
    #println("length",length(setdiff(links, core2nodes)))
    #println("current size:", nbc2)
    println("Number of iterations: ", length(links))
    for (idx,l) in enumerate(links)
        if tms[linksneighbor[idx]] > 1
            push!(bridges,l)
            push!(bridges,linksneighbor[idx])
        else
            #println(idx," ( ", l," ", linksneighbor[idx], " )")
            p,d = bfs_edge_subtree2(g, l, linksneighbor[idx], tms)        
            if length(d) == 1
                ## TODO: assert length(path) == 2
                push!(bridges,p[1])
                push!(bridges,p[2])
            else 
                cntr = zeros(Int64, maximum(d))
                cntr = count2(d,length(d))
                if isempty(findin(core2nodes, l))
                    #println(l, " ", size(cntr,1))
                    #sizes[nbc2 + i] = [];
                    for j in 1:length(cntr)
                        #println("idx:", nbc2 + i, " ->", cntr[j])
                        #println(sizes[nbc2 + i])
                        #push!(sizes[nbc2 + i], cntr[j])
                        #println(sizes[nbc2 + i])
                    end
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

#println("final bridges = ", bridges)
#println("core2nodes:", core2nodes)
#println("sizes:", sizes)
#for i in 1:length(sizes)
#    println(sizes[i])
#end

B = Bridges(bridges, core2nodes, sizes)
println("### cleaning up bridges time:",time()- t, "(s)")
t = time()
shell1nodes = Array{Int64,1}
shell1nodes = k_shell(g, 1)
#println(length(shell1nodes))
println("### finding shell 1 time:",time()- t, "(s)")
return B, shell1nodes
end

### similar to bridges2 but with different return type.
### It is used in approxcore2 method and it only returns
### the bdry nodes and their related paths of core1.
### It doesn't keep track of the actual bridges between
### components of size > 1.
function bridges3(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
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
    #println(sizes)
    

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
                sizes[nbc2 + i] = [];
                for j in 1:length(cntr)
                    push!(sizes[nbc2 + i], cntr[j])
                end
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

#println("core2nodes:", length(core2nodes))
#println("sizes:", size(sizes,1))
# for i in 1:length(sizes)
#     println(sizes[i])
# end
# println("final bridges = ", bridges)

println("COUNT TIMING IS ",time()- t, "(s)")

return core2nodes, sizes
end


### TODO start BFS from linksneighbor
### TODO stop if node has tms > 1 or if length(outneighbors(g, v) == 1
function count2( array ::  Array{Int64,1}, len :: Int64)
    cntr = zeros(Int64, maximum(array))
    #cntr =  Array{Int64,1}()
    #println(array)
    c :: Int64 = 1
    i :: Int64 = 1
    for (idx,u) in enumerate(array)
        if idx > 1
            if u == array[idx - 1]
                c += 1
            else
                cntr[i] = c
                c = 1
                i += 1
            end
        end
        if idx == len
            cntr[i] = c

        end
    end
    return cntr
end

### same as below, but here the size of the array
### that holds the returned values is maximum(array)
function count3( array ::  Array{Int64,1}, len :: Int64)
    cntr = zeros(Int64, maximum(array))
    #cntr =  Array{Int64,1}()
    #println(array)
    c :: Int64 = 1
    i :: Int64 = 1
    for (idx,u) in enumerate(array)
        if idx > 1
            if u == array[idx - 1]
                c += 1
            else
                cntr[i] = c
     
                c = 1
                i += 1
            end
        end
        if idx == len
            cntr[i] = c
     
        end
    end
    return cntr
end



### count how many times a value appears in an array
### the size of the array to store this info is equal to
### unique(array)
function count4( array ::  Array{Int64,1}, len :: Int64)
    cntr =  Array{Int64,1}()
    c :: Int64 = 1
    for (idx,u) in enumerate(array)
        if idx > 1
            if u == array[idx - 1]
                c += 1
            else
                push!(cntr, c)
                c = 1
            end
        end
        if idx == len
            push!(cntr, c)
     
        end
    end
    return cntr
end




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



### following function is inspired by bfs_edge_subtree of lightGraphs
### but has slightly different functionality.
### function to start bfs from edge (source,next)
function bfs_edge_subtree3(g::AbstractGraph{T}, source :: Int64, next:: Int64, tms ::  Array{Int64,1} ) where T

    #println(next)
    if tms[next] > 1
        return [next], [1], []
    end
    n = nv(g)
    parents = zeros(T, nv(g))
    visited = falses(n)
    distance = Vector{T}()
    path = Vector{T}()
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    visited[source] = true
    visited[next] = true
    #parents[source] = source
    push!(cur_level, next)
    cnt :: Int64 = 1
    #push!(distance,0)
    push!(distance,cnt)
    push!(path,next)

    cnt += 1
    extra = []
   
    while !isempty(cur_level)
        @inbounds for v in cur_level
# check what @simd is for
            @inbounds @simd for i in  outneighbors(g, v)
                if !visited[i] && tms[i] == 1
                    push!(next_level, i)
                    push!(distance,cnt)
                    push!(path,i)
                    parents[i] = v
                    visited[i] = true
                elseif !visited[i] && tms[i] > 1
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
                    # I do not push i to the next_level
                    # so none of its neighbors can be reached
                    push!(path,i)
                    parents[i] = v
                    visited[i] = true
                    push!(extra,i)
                end
            end
        end
        empty!(cur_level)
        cnt += 1
        cur_level, next_level = next_level, cur_level
        ## do I need sorting?
        #sort!(cur_level)
    end
    if isempty(extra) == true
        return path, distance, extra
    else
        return path, distance, parents
    end
end




### function to start bfs from edge (source,next)
function bfs_edge_subtree2(g::AbstractGraph{T}, source :: Int64, next:: Int64) where T
    
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
        ## do I need sorting?
        #sort!(cur_level)
    end
    return path, distance
end

####################################################################

function locateBridges(A :: SparseMatrixCSC{Float64})
    edges = bridges(LightGraphs.Graph(A))
    A, extnodes = removeBridges(A, edges, nedges)
    B = Bridges(edges, extnodes)
    return A,B
end

function removeBridges(A :: SparseMatrixCSC{Float64}, brs, nbrs :: Integer)
    nodes =  Set{Int64}()
    for e in brs
        A[e.src,e.dst] = 0.0
        A[e.dst,e.src] = 0.0
        push!(nodes,e.src)
        push!(nodes,e.dst)
    end
    return dropzeros!(A), nodes
end

function removeBridges(A :: SparseMatrixCSC{Float64}, B :: Bridges, core1nodes :: Array{Int64,1})
    ### delnodes creates a new array so the numbering
    ### is remapped. I want the numbering to stay the
    ### same, so I will just write zeros ontop of
    ### A whereever is needed.
    ### A = delnodes(A, B.core1nodes)
    j :: Int64 = 0
    rows = rowvals(A)
    ## or   colptr = mat.colptr , rowval = mat.rowval
    #vals = nonzeros(A)
    for u in core1nodes
        j = nzrange(A, u)[1]
        A[u,rows[j]] = 0.0
        A[rows[j],u] = 0.0
    end
    for i in 1:2:length(B.edges)
        A[B.edges[i],B.edges[i+1]] = 0.0
        A[B.edges[i+1],B.edges[i]] = 0.0
    end
    dropzeros!(A)
    #println("size of A after dropzeros! = ",size(A,1)," ",nnz(A))
    return A
end

function removeBridges(A :: SparseMatrixCSC{Float64}, B :: Bridges, core1nodes :: Set{Int64})
    j = Int64
    rows = rowvals(A)
    for u in core1nodes
        j = nzrange(A, u)[1]
        A[u,rows[j]] = 0.0
        A[rows[j],u] = 0.0
    end
    for i in 1:2:length(B.edges)
        A[B.edges[i],B.edges[i+1]] = 0.0
        A[B.edges[i+1],B.edges[i]] = 0.0
    end
    return dropzeros!(A)
end

function extractBridges(A :: SparseMatrixCSC{Float64})

    B = Bridges
    #start_time = time()
    #B, core1nodes = bridges5(LightGraphs.Graph(A))
    #println((100*B.m)/(nnz(A)/2), "% edges are bridges type core2.")
    #println(100*length(B.edges)/A.n, "% nodes are core2.")
    #println("### finding bridges time (bridges5): ", time() - start_time, "(s)")
    #start_time = time()
    #B, core1nodes = bridges6(LightGraphs.Graph(A))
    # #println((100*B.m)/(nnz(A)/2), "% edges are bridges type core2.")
    # #println(100*length(B.edges)/A.n, "% nodes are core2.")
    #println("### finding bridges time (bridges6): ", time() - start_time, "(s)")
    println()
    start_time = time()
    B, core1nodes = bridges6(LightGraphs.Graph(A))
    # #println((100*B.m)/(nnz(A)/2), "% edges are bridges type core2.")
    # #println(100*length(B.edges)/A.n, "% nodes are core2.")
    println("### finding bridges time (bridges6): ", time() - start_time, "(s)")
    println()
    # println("### Time for bridges2 ")
    start_time = time()
    B, core1nodes = bridges2(LightGraphs.Graph(A))
     println("### finding bridges time (bridges2): ", time() - start_time, "(s)")
    # # println("### Time for bridges5 ")
    t = time()
    A  = removeBridges(A, B, core1nodes)
    println("### remove bridges time: ", time()- t, "(s)")
    return A,B
end


function computeCore2Bridges(A :: SparseMatrixCSC{Float64})
    g = LightGraphs.Graph(A)
    Bridgescore1 :: Int64 = 0
    Bridgescore2 :: Int64 = 0
    Bridgescore1 = nbofbridges(g)
    v1 = k_core(g,2)
    A1 = A[v1,v1]
    Bridgescore2 = nbofbridges(LightGraphs.Graph(A1))
    println("Bridgescore1 = ", Bridgescore1, "Bridgescore2 = ", Bridgescore2)
    println((100*Bridgescore2)/(nnz(A)/2), "% edges are bridges type core2.")
end




### testing bridges function that is used to calculate bridges
### the bndry nodes and all the related to bndry nodes
### simple paths of core1 type. 
function bridges4(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    m :: Int64 = 0
    tms = zeros(T, nv(g))
    degree1neighbor = Array{Int64,1}()

    st = time()
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
    
    println("### finding all bridges time:",time()- st, "(s)")

    t = time()

    uedges = unique(edges)
    m = length(edges)
    
    candidates = Array{Int64,1}()
    extra = Array{Int64,1}()
    shellone = Array{Int64,1}()
    ### create initial groups of nodes, candidates for core2nodes,
    ### extra nodes to be tested and shell-one nodes
    for i in uedges
        if tms[i] > 1
            push!(extra,i)
        end
        tms[i] = length(outneighbors(g, i)) - tms[i] + 1
        if tms[i] > 1
            push!(candidates,i)
        else
            push!(shellone,i)
        end
    end
    ### create the disconnected subgraph
    ### TODO: I know the length! m/2
    elist = Edge{Int64}[]
    for i in 1:2:m
        v = edges[i]
        w = edges[i+1]
        edge = v < w ? Edge(v, w) : Edge(w, v)
        push!(elist, edge)
    end
    sg, vmap = induced_subgraph(g, elist)
    #println(sg)
    #println(vmap)
    bridges =  Array{Int64,1}()
    nbCandidates = length(candidates)
    #println("Candidates:")
    #println(candidates)
    
    sizes = Array{Array{Int, 1}}(nbCandidates)
    for i in indices(sizes,1) sizes[i] = []; end
    visited = falses(nbCandidates) ## same size as vmap
    removeCandidate =  Array{Int64,1}()
    
    for (i,u) in enumerate(candidates)
        umap = getindex(findin(vmap,u))
        family , degree = bfs_node_subtree2(sg, umap)
        #println("node $u has path", family)
        #println("node $u has dist", degree)
        next = family[findin(degree,1)]
        ### remove from candidates of core2nodes
        ### all vertices that belong to bridges (second condition)
        ### and all vertices that belong to core 1 (first condition)
        if tms[u] == 1 || (tms[u] > 1 && tms[vmap[next[1]]] > 1)
            push!(removeCandidate,i)
        end
        for (idx,v) in enumerate(next)
            if tms[vmap[v]] > 1
                if visited[i] == false
                    push!(bridges,u)
                    push!(bridges,vmap[v])
                    visited[getindex(findin(candidates,u))] = true
                    visited[getindex(findin(candidates,vmap[v]))] = true
                end
                #d[idx] = 0
                break
            end
        end
        ### calculate sizes
        sizes[i] = count3(degree,length(degree))
        #println("node $u has sizes", sizes[i])
    end
    #println("Dealing with extra nodes!!!!")
    extra = setdiff(extra, candidates)
    #println(extra)
    evisited = falses(length(extra)) ## same size as extra
    sizes2 = Array{Array{Int, 1}}(length(extra))
    for (i,u) in enumerate(extra)
        umap = getindex(findin(vmap,u))
        evisited[i] == true;
        family , degree = bfs_node_subtree2(sg, umap)
        #println("node $u has path", family)
        #println("node $u has dist", degree)
        cnt :: Int64 = 0
        ngbors = Array{Int64, 1}()
        for (idx,v) in enumerate(family)
            if tms[vmap[v]] > 1
                cnt += 1
                push!(ngbors,v)
            end
        end
        if cnt <= 1
            continue;
        else
            ## TODO: remove ngbors form  shellone
            removeShellone = Array{Int64,1}()
            #println("node $u should be in a new component!!")
            #println("with bridges to: ",vmap[ngbors], ". We don't know if direct edges yet!")
            for (j,v) in enumerate(ngbors)
                if has_edge(g, u, vmap[v]) == true
                    push!(bridges,vmap[v])
                    push!(bridges,u)
                    push!(removeShellone,u)
                    push!(removeShellone,vmap[v])
                    ### I have to correct sizes of connecting nodes:
                    
                else
                    println("WARNING: have to deal with path of nodes (not simple edge)! ")
                    path, distance = bfs_node_subtree2(sg, umap , v)
                    # for p in path[1:end-1]
                    #     visited[getindex(findin(extra,vmap[p]))] = true
                    # end
                    # next = p[1]
                    # push!(bridges,vmap[next])
                    # push!(bridges,u)
                    # push!(removeShellone,u)
                    # push!(removeShellone,vmap[p])
                    #### TODO: uncomment the above
                    
                end
                # println(path[length(path)])
                # println(vmap[path[length(path)]])
                # push!(bridges,vmap[path[length(path)]])
                # push!(bridges,u)

                #sizes2[i] = count3(degree,length(degree))
                #println("node $u has sizes", sizes[i])
            end
            idx = findin(shellone,removeShellone)
            #println("deleting shellones: ",removeShellone, idx)
            deleteat!(shellone,idx)
            
        end
    end
    #println("Dealing with extra nodes- Over!!!!")
    
    #println(removeCandidate)
    deleteat!(candidates,removeCandidate)
    deleteat!(sizes,removeCandidate)
    nbc2 :: Int64 = length(candidates) 
    # for (i,u) in enumerate(candidates)
    #     println(sizes[i])
    # end
    #println("Final:")
    #println("bridges:", bridges) ## bridges
    #println("core2nodes",candidates)  ## core2nodes
    #println("sizes",sizes)  ## core2nodes
    #println("shellone",shellone)

#println("core2nodes:", core2nodes)
#println("sizes:", sizes)
# for i in 1:length(sizes)
#     println(sizes[i])
# end
    # println("final bridges = ", bridges)
core2nodes = Array{Int64,1}(nbc2)
core2nodes = candidates
if (length(sizes) != length(core2nodes))
    println("WARNING: Wrong size of sizes and core2nodes!")
    println("WARNING: Program will exit!")
end
B = Bridges(bridges, core2nodes, sizes)
return B, shellone
end

function bfs_node_subtree2(g::AbstractGraph{T}, source) where T    
    n = nv(g)
    visited = falses(n)
    distance = Vector{T}()
    path = Vector{T}()
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    visited[source] = true
    push!(cur_level, source)
    cnt :: Int64 = 1
    
    push!(path,source)
    
    while !isempty(cur_level)
        @inbounds for v in cur_level
# check what @simd is for
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
        ## do I need sorting?
        #sort!(cur_level)
    end
    return path, distance
end



function bfs_node_subtree2(g::AbstractGraph{T}, source , dest) where T    
    n = nv(g)
    visited = falses(n)
    distance = Vector{T}()
    path = Vector{T}()
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    visited[source] = true
    push!(cur_level, source)
    cnt :: Int64 = 1
    #println("starting from: ", source)
    #println("destination: ", dest)
    #push!(path,source)
    while !isempty(cur_level)
        @inbounds for v in cur_level
# check what @simd is for
            @inbounds @simd for i in  outneighbors(g, v)
                if i == dest
                    #push!(path,i)
                    #println("Final path:", path)
                    return path, distance
                end
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
        #println(path)
        ## do I need sorting?
        #sort!(cur_level)
    end
    return path, distance
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
                    println("WARNING: code may has some logical bugs!")
                    push!(distance,cnt)
                    last = path[end]
                    println("WARNING: the bridges that are not treated as bridges:",path)
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




function verifyBridge(next, tms, visited, vmap, bridges, candidates, u ,i)
    for (idx,v) in enumerate(next)
        if tms[vmap[v]] > 1
            if visited[i] == false
                push!(bridges,u)
                push!(bridges,vmap[v])
                visited[getindex(findin(candidates,u))] = true
                visited[getindex(findin(candidates,vmap[v]))] = true
            end
            #d[idx] = 0
            break
        end
    end
end


### testing bridges function that is used to calculate bridges
### the bndry nodes and all the related to bndry nodes
### simple paths of core1 type. 
function bridges5(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    m :: Int64 = 0
    tms = zeros(T, nv(g))
    degree1neighbor = Array{Int64,1}()

    st = time()
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
    
    println("### finding all bridges time:",time()- st, "(s)")
    
    t = time()

    uedges = unique(edges)
    m = length(edges)
    candidates = Array{Int64,1}()
    stones = Array{Int64,1}()
    peripheral = Array{Int64,1}()
    ### create groups of nodes, candidates for core2nodes,
    ### stones for intermediate nodes and peripheral
    ### for degree one nodes.
    for i in uedges
        tms[i] = length(outneighbors(g, i)) - tms[i] + 1
        if tms[i] > 1
            push!(candidates,i)
        else
            if length(outneighbors(g, i)) == 1
                push!(peripheral,i)
            else
                push!(stones)
            end
        end
    end
    
    shellone =  [peripheral; stones]
    #### TODO: are stones always empty??
    if isempty(stones) == true
        println("WARNING: stones are empty!!")
    end
    #println("stones: ", stones)
    #println("candidates:", candidates)
    #println("peripheral:", peripheral)
    
    sg, vmap = createForest(g, edges, m)
    localPeripheral = Array{Int64, 1}(length(peripheral))
    localPeripheral = findin(vmap,peripheral)
    ### maybe detectCore2 should return a reduced graph!
    ### because sg has still a large number of vertices.
    ### (In detectCore2 edges are removed instead of vertices
    ### to avoid renumbering vertices for each iteration.)
    ### This may have a bad impact in the performance
    ### of detectBridges
    #sg, sizes= detectCore(sg, localPeripheral, vmap, candidates)
    sg, sizes= detectCore2(sg, localPeripheral, vmap, candidates)
    core2nodes = Array{Int64,1}()
    bridges =  Array{Int64,1}()
    ### detectBridges finds the obvious bridges, that
    ### could be found directly from uedges, but
    ### it also finds those extra bridges that are needed
    ### if a path of nodes connects two components.
    bridges, core2nodes= detectBridges(sg, vmap,sizes)
    
    println("Final:")
    println("bridges:", vmap[bridges]) ## bridges
    println("core2nodes",vmap[core2nodes])  ##core2nodes
    for i in sizes
         if !isempty(i)
             print(i, " ")
         end
     end
     println()
    #println("shellone", shellone)
    uedges = setdiff(uedges,vmap[bridges])
    uedges = setdiff(uedges,vmap[core2nodes])
    #println("shellone", uedges)
    #println("len sizes =", size(sizes,1), " len core2nodes =", length(core2nodes))
    # B = Bridges(vmap[bridges], vmap[core2nodes], sizes)
    new_sizes = Array{Array{Int64,1}}(length(core2nodes))
    j :: Int64 = 1
    for i in sizes
        if !isempty(i)
            new_sizes[j] = i
            j += 1
        end
    end
if length(new_sizes) != length(core2nodes)
    println("WARNING: Wrong size of sizes and core2nodes!")
    println("WARNING: Program will exit!")
    exit()
end

if length(uedges) != (length(core2nodes) + length(bridges)) 
    println("WARNING: something's wrong with bridges!")
end

B = Bridges(vmap[bridges], vmap[core2nodes], new_sizes)
return B, uedges
end


function createForest(g::AG, edges:: Array{Int64,1}, m:: Int64)  where {T, AG<:AbstractGraph{T}}
    ### create the disconnected subgraph
    ### TODO: I know the length! m/2

    if isodd(m) == true
        println("WARNING: the number of nodes adjacent to bridges should be even!")
        println("WARNING: program will exit!")
    end
    elist = Edge{Int64}[]
    #elist = Edge{Int64}(m/2)
    for i in 1:2:m
        v = edges[i]
        w = edges[i+1]
        edge = v < w ? Edge(v, w) : Edge(w, v)
        push!(elist, edge)
    end
    sg, vmap = induced_subgraph(g, elist)
    return sg, vmap
end


function detectCore(sg::AG, localPeripheral :: Array{Int64,1}, vmap:: Array{Int64,1},candidates:: Array{Int64,1})  where {T, AG<:AbstractGraph{T}}
    
    sn = length(vmap)
    #println(sn)
    visited = falses(sn) ## same size as vmap
    visited[findin(vmap, candidates)] = true
    # println("visited",visited)
    sizes = Array{Array{Int, 1}}(sn)
    for i in indices(sizes,1) sizes[i] = []; end
    j :: Int64 = 1
    while !isempty(localPeripheral)
        #println(sg)
        if length(localPeripheral) < 20
            #println("Peripheral",vmap[localPeripheral])
            #for u in localPeripheral
              #  println("u ",vmap[u], " next ", vmap[outneighbors(sg, u)])
             #end
        end
        count = zeros(Int64,sn)
        for (idx, node) in enumerate(localPeripheral)
            next = outneighbors(sg, node)[1]
            #if length(localPeripheral) < 20
            #    println("node ",vmap[node], " next ", vmap[next])
            #end
            count[next] += 1
            ### TODO: this part, not so correct!
            if isempty(sizes[node]) == false
                for (idx, s) in enumerate(sizes[node])
                    #if isempty(sizes[next][idx]) == false
                    #    println(sizes[next][idx])
                    #    #sizes[next][idx] += s
                    #else
                        push!(sizes[next], s)
                    #end
                end                
                # for i in 1:length(sizes[next])
                #     println("$i , sizes:",sizes[i])
                # end

                #if length(localPeripheral) < 20
                    #println("next ",vmap[next], " is set to: ", sizes[next])
                    #println("node ",vmap[node], " before ", sizes[node])
                #end
                #pop!(sizes[node])
                empty!(sizes[node])
                #if length(localPeripheral) < 20
                    #println("node ",vmap[node], " removed ", sizes[node])
                #end
            end
            rem_edge!(sg, node, next)
            
            visited[node] = true
        end
        for (i,value) in enumerate(count)
            if value != 0
                push!(sizes[i], value)
            end
        end
        # println("New level size:")
        # for (i,s) in enumerate(sizes)
        #     if !isempty(s)
        #         l = vmap[i]
        #         println("node $l = ", s)
        #     end
        # end

        # println(visited)
        empty!(localPeripheral)
        for u in vertices(sg)
            if !visited[u] && length(outneighbors(sg, u)) == 1
                push!(localPeripheral,u)
            end
        end
        j +=1
    end
    return sg, sizes
end


function detectCore2(sg::AG, localPeripheral :: Array{Int64,1}, vmap:: Array{Int64,1},candidates:: Array{Int64,1})  where {T, AG<:AbstractGraph{T}}
    
    sn = length(vmap)
    visited = falses(sn) ## same size as vmap
    visited[findin(vmap, candidates)] = true
    sizes = Array{Array{Int, 1}}(sn)
    for i in indices(sizes,1) sizes[i] = []; end
    while !isempty(localPeripheral)
        #println(sg)
        #if length(localPeripheral) < 20
        #println("Peripheral",vmap[localPeripheral])
        #for u in localPeripheral
        #  println("u ",vmap[u], " next ", vmap[outneighbors(sg, u)])
        #end
        #end
        count = zeros(Int64,sn)
        for (idx, node) in enumerate(localPeripheral)
            next = outneighbors(sg, node)[1]
            count[next] += 1
            ### make room for the edge connection
            if isempty(sizes[next]) == true
                push!(sizes[next], 0)
            end
            ### update path connections
            len = length(sizes[next])
            for (i, s2) in enumerate(sizes[node])
                if i+1 <= len
                    sizes[next][i+1] += s2
                else
                    push!(sizes[next],s2)
                end
            end
           
            empty!(sizes[node])
            rem_edge!(sg, node, next)           
            visited[node] = true
        end
        ### adding the edge connection, alway in the first place
        for (i,value) in enumerate(count)
            if value != 0
                if length(sizes[i]) > 0
                    sizes[i][1] += value
                else
                    #### TODO: remove this if else.
                    #### length of sizes[i] when value
                    #### is not zero is always > 0 
                    println("NEVER HERE!")
                    push!(sizes[i], value)
                end
            end
        end        
        # println("New level size:")
        # for (i,s) in enumerate(sizes)
        #     if !isempty(s)
        #         l = vmap[i]
        #         println("node $l = ", s)
        #     end
        # end
        empty!(localPeripheral)
        for u in vertices(sg)
            if !visited[u] && length(outneighbors(sg, u)) == 1
                push!(localPeripheral,u)
            end
        end
    end
    return sg, sizes
end


function detectBridges(sg::AG, vmap:: Array{Int64,1},sizes :: Array{Array{Int, 1}})  where {T, AG<:AbstractGraph{T}}
    core2nodes = Array{Int64,1}()
    bridges =  Array{Int64,1}()
    sn = length(vmap)
    visited = falses(sn)
    
    for (idx, u) in enumerate(vertices(sg))
        if !isempty(sizes[idx])
            push!(core2nodes,u)
        end
        if isempty(outneighbors(sg, u))
            continue; 
        elseif visited[u] == false
            next, dist = bfs_node_subtree2(sg, u)
            #println("distance:", dist)
            ### no branching!
            if dist[end] == length(next) - 1
                addBridgesInPath!(next,bridges,visited,vmap,dist)
                #### branching!
            elseif dist[end] < length(next) - 1
                addBridgesInTree!(next,bridges,visited,vmap,dist)
            end
        end
    end
    return bridges, core2nodes
end


function detectBridges2(sg::AG, vmap:: Array{Int64,1},sizes :: Array{Array{Int, 1}})  where {T, AG<:AbstractGraph{T}}
    core2nodes = Array{Int64,1}()
    bridges =  Array{Int64,1}()
    sn = length(vmap)
    visited = falses(sn)
    
    for (idx, u) in enumerate(vertices(sg))
        if !isempty(sizes[idx])
            push!(core2nodes,u)
        end
        if isempty(outneighbors(sg, u))
            continue; 
        elseif visited[u] == false
            next, dist = bfs_node_subtree2(sg, u)
            #println("distance:", dist)
            ### no branching!
            if dist[end] == length(next) - 1
                addBridgesInPath!(next,bridges,visited,vmap,dist)
                #### branching!
            elseif dist[end] < length(next) - 1
                addBridgesInTree2!(next,bridges,visited,vmap,dist)
            end
        end
    end
    return bridges, core2nodes
end

function addBridgesInPath!(next:: Array{Int64,1}, bridges:: Array{Int64,1}, visited :: BitArray{1},vmap :: Array{Int, 1}, distances :: Array{Int64,1})
    #println("WARNING: pathssssss!! But no branching!!")    
    #println("u ", vmap[next[1]], " -->", vmap[next], "--> dist", distances)
    if size(next,1) == 2
       
        push!(bridges, next[1])
        push!(bridges, next[2])
        visited[next[1]] = true
        visited[next[2]] = true
    elseif size(next,1) > 2
        push!(bridges, next[1])
        push!(bridges, next[2])
        visited[next[1]] = true
        visited[next[2]] = true
        push!(bridges, next[end])
        push!(bridges, next[end-1])
        visited[next[end]] = true
        visited[next[end-1]] = true               
    end
    return distances, visited
end


function addBridgesInTree!(next:: Array{Int64,1}, bridges:: Array{Int64,1}, visited :: BitArray{1},vmap :: Array{Int, 1}, distances :: Array{Int64,1})
    println("WARNING: pathssssss with Branching!!")    
    m = maximum(distances)
    #println("u ", vmap[next[1]], " -->", vmap[next])
    #println("u ", vmap[next[1]]," --> dist", distances)
    #println(count3( distances, length(distances)))
    #println(count2( distances, length(distances)))
    count = count3( distances, length(distances))
    k = 0
    for (i, c) in enumerate(count)
        if c > 1
            for j in 1:c
                #println(vmap[next[k+1+j]])
                push!(bridges, next[k+i+j])
                visited[next[k+i+j]] = true
                push!(bridges, next[1])
                visited[next[1]] = true
            end
            break;
        end
        k += c
    end
    # push!(bridges, next[1])
    # push!(bridges, next[2])
    # visited[next[1]] = true
    # visited[next[2]] = true
    # for (idx, u) in enumerate(next)
    #     if distances[idx] == m
    #         push!(bridges, u)
    #         push!(bridges, next[idx-1])
    #         visited[next[u]] = true
    #         visited[next[idx-1]] = true               
    #     end
    # end
    #TODO: uncomment this!!
    return distances, visited
end

### add the end points not the first
function addBridgesInTree2!(next:: Array{Int64,1}, bridges:: Array{Int64,1}, visited :: BitArray{1},vmap :: Array{Int, 1}, distances :: Array{Int64,1})
    println("WARNING: pathssssss with Branching!!")    
    m = maximum(distances)
    println("u ", vmap[next[1]], " -->", vmap[next])
    println("u ", vmap[next[1]]," --> dist", distances)
    println(count3( distances, length(distances)))
    #println(count2( distances, length(distances)))
    count = count3( distances, length(distances))
    k = 0
    for (i, c) in enumerate(count)
        if i == length(count)
            for j in 1:c
                println(vmap[next[k+1+j]])
                push!(bridges, next[k+1+j])
                visited[next[k+1+j]] = true
                ##### TODO: add the other node of the bridge #####
                ##### which is in the previous level ###### 
            end
            break;
        end
        k += c
    end
    return distances, visited
end


### main bridges function that is used to calculate bridges
### the bndry nodes and all the related to bndry nodes
### simple paths of core1 type. 
function bridges6(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    m :: Int64 = 0
    tms = zeros(T, nv(g))
    degree1 = Array{Int64,1}()
    st = time()
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
                    push!(edges,w)
                    tms[w] += 1
                    push!(edges,v)
                    tms[v] += 1
                    if length(v_neighbors) == 1
                        push!(degree1, w)
                    elseif length(outneighbors(g, w)) == 1
                        push!(degree1, v)
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
    
    println("### finding all bridges time:",time()- st, "(s)")
    t = time()
    degree1 = sort(degree1)
    udegree1count = count4(degree1, length(degree1))
    udegree1 = unique(degree1)

    core1nodes = Array{Int64,1}()
    core2nodes = Array{Int64,1}()
    uedges = unique(edges)
    
    for i in uedges
        tms[i] = length(outneighbors(g, i)) - tms[i] + 1
        if tms[i] > 1
            push!(core2nodes, i)
        else
            push!(core1nodes, i)
        end
    end
    
    remove = Array{Int64,1}()
    
    for (idx, u) in enumerate(udegree1)
        if tms[u] == 1
            push!(remove,idx)
        end
    end
    deleteat!(udegree1,remove)
    deleteat!(udegree1count,remove)
    
    if core2nodes == unique(core2nodes)
        println("core2nodes unique!!")
    else
        println("WARNING!!")
        exit()
    end
    if udegree1 == unique(udegree1)
        println("udegree1 unique!!")
    else
        println("WARNING!!")
        exit()
    end
    if core1nodes == unique(core1nodes)
        println("core1nodes unique!!")
    else
        println("WARNING!!")
        exit()
    end
    
    links = Array{Int64,1}()
    linksneighbor = Array{Int64,1}()
    
    nbc2 :: Int64 = length(core2nodes)
    position = findin(edges,core2nodes)
    
    ### links are not unique but linksneighbor are unique!!
    ### CHECK: the above statement to make sure!
    ### links, linksneighbor correspond to edges of core2node
    ### links and linksneighbor have the same size, but
    ### core2node is basically unique(links) and should have
    ### the same size as the first dimension of sizes.  
    for (idx, i) in enumerate(position)
        if isodd(i) == true
            if tms[edges[i]] > 1 && length(outneighbors(g, edges[i+1])) > 1
                push!(links, edges[i])
                push!(linksneighbor, edges[i + 1])
            end
        else
            if tms[edges[i]] > 1 && length(outneighbors(g, edges[i-1])) > 1
                push!(links, edges[i])
                push!(linksneighbor, edges[i - 1])
            end
        end
    end
    
    ### at this point we have:
    ### in core1nodes all nodes with tms 1.
    ### It means that we should remove some of them (those that are present in bridges).
    ### In core2nodes we have the opposite:
    ### we have all nodes that belong to a bigger component
    ### but we may need to add core2nodes that are part
    ### of smaller component (or even single-node component)
    ### but we don't know it yet
    
    sizes = Array{Array{Int, 1}}(nbc2)
    for i in indices(sizes,1) sizes[i] = []; end

    bridges =  Array{Int64,1}()
#### TODO: maybe visited should be of size all edges of bridges type
#### to avoid duplicate calls to bfs_edge_subtree3!!
    visited = falses(length(links))
    boundary = Array{Int64,1}()

    println("Number of iterations: ", length(links))
    println("smaller number of iter ? : ", length(core2nodes))
    println("### setting up time:",time()- t, "(s)")    
    tr1 = time()
    for (i, u) in enumerate(links)
        v = linksneighbor[i]    
        ### this is directly a bridge
        if tms[v] > 1
            ## avoid double entries in bridges
            if (getindex(visited[i]) == false &&
                getindex(visited[findin(links, v)]) == false)
                push!(bridges, u)
                visited[i] = true
                push!(bridges, v)
                visited[findin(links, v)] = true
            end
        else
            tree, distances, extra = bfs_edge_subtree3( g, u, v, tms)
            if isempty(extra) == true
                unii = getindex(findin(core2nodes, u))
                len = length(sizes[unii])
                cntr = zeros(Int64, maximum(distances))
                cntr = count3(distances,length(distances))
                if len == 0
                    push!(boundary, u)
                end
                for (j,c) in enumerate(cntr)
                    if j <= len 
                        sizes[unii][j] += c
                    else
                        push!(sizes[unii], c)
                    end 
                end
            else
                #### means that we have intermediate nodes to
                #### be single node components!
                #### we make them all singular components
                #### TODO : I need a visited here!!!
                for i in 1:length(extra)
                    if extra[i]!=0
                        push!(bridges, i)
                        push!(bridges, extra[i])
                    end
                end
                push!(bridges, u)
                push!(bridges, tree[1])
            end
        end
    end
println("### exploration time:",time()- tr1, "(s)")

tr3 = time()
### make sure that bridges are not part of core1nodes
deleteat!(core1nodes,findin(core1nodes, bridges))

if boundary == unique(boundary)
    println("boundary unique!!")
else
    println("WARNING!!")
    exit()
end


additional = setdiff(udegree1, boundary)

new_sizes = Array{Array{Int64,1}}(length(boundary) + length(additional))
for i in indices(new_sizes,1) new_sizes[i] = []; end


j :: Int64 = 1
for i in 1:length(sizes)
    if !isempty(sizes[i])
        new_sizes[j] = sizes[i]
        u = core2nodes[i]
        if in(u, udegree1) == true 
            new_sizes[j][1] += getindex(udegree1count[findin(udegree1, u)])
        end
        j += 1
    end
end

### updating additional nodes, that only have degree1 connections
for (i, u) in enumerate(udegree1)
    if in(u, additional) == true
        new_sizes[j] = [udegree1count[i]]
        j += 1
    end
end

boundary = [boundary ; additional]

println("bridges = ",bridges)
println("bdry = ", boundary)    
println("sizes = ",new_sizes)
#println("core1nodes = ",core1nodes)



if length(new_sizes) != length(boundary)
    println("WARNING: Wrong size of sizes and core2nodes!")
    println("WARNING: Program will exit!")
    exit()
end
for u in core1nodes
    if in(u,boundary) == true || in(u,bridges) == true
        println("WARNING: core1nodes are not correct!")
        exit()
    end
end

println("### reshaping bridges time:",time()- tr3, "(s)")

B = Bridges(bridges, boundary, new_sizes)
println("### cleaning up time:",time()- t, "(s)")
return B, core1nodes
end


function getCoreNodes(g::AG, uedges :: Array{T,1}, tms :: Array{T,1}) where {T, AG<:AbstractGraph{T}}
    core1nodes = Array{T,1}()
    core2nodes = Array{T,1}()
    for i in uedges
        tms[i] = length(outneighbors(g, i)) - tms[i] + 1
        if tms[i] > 1
            push!(core2nodes, i)
        else
            push!(core1nodes, i)
        end
    end
    return core1nodes, core2nodes
end

function cleanDegree1Nodes!(g::AG, tms :: Array{T,1}, udegree1 :: Array{T,1}, udegree1count :: Array{T,1}) where {T, AG<:AbstractGraph{T}}
    
    remove = Array{Int64,1}()
    for (idx, u) in enumerate(udegree1)
        if tms[u] == 1
            push!(remove,idx)
        end
    end
    deleteat!(udegree1,remove)
    deleteat!(udegree1count,remove)
    return udegree1,udegree1count
end


### main bridges function that is used to calculate bridges
### the bndry nodes and all the related to bndry nodes
### simple paths of core1 type. 
function bridges7(g::AG) where {T, AG<:AbstractGraph{T}}
    s = Vector{Tuple{T, T, T}}()
    low = zeros(T, nv(g)) #keeps track of the earliest accesible time of a vertex in DFS-stack, effect of having back-edges is considered here
    pre = zeros(T, nv(g)) #checks the entry time of a vertex in the DFS-stack, pre[u] = 0 if a vertex isn't visited; non-zero, otherwise
    # bridges = Edge{T}[]   #keeps record of the bridge-edges
    edges = Array{Int64,1}()
    m :: Int64 = 0
    tms = zeros(T, nv(g))
    degree1 = Array{Int64,1}()
    st = time()
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
                    push!(edges,w)
                    tms[w] += 1
                    push!(edges,v)
                    tms[v] += 1
                    if length(v_neighbors) == 1
                        push!(degree1, w)
                    elseif length(outneighbors(g, w)) == 1
                        push!(degree1, v)
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
    
    println("### finding all bridges time:",time()- st, "(s)")
    t = time()
    core1nodes = Array{Int64,1}()
    core2nodes = Array{Int64,1}()
    degree1 = sort(degree1)
    udegree1count = count4(degree1, length(degree1))
    udegree1 = unique(degree1)
    
    remove = Array{Int64,1}()
    uedges = unique(edges)
    m = length(edges)
    
    for i in uedges
        tms[i] = length(outneighbors(g, i)) - tms[i] + 1
        if tms[i] > 1
            push!(core2nodes, i)
        else
            push!(core1nodes, i)
        end
    end

    for (idx, u) in enumerate(udegree1)
        if tms[u] == 1
            push!(remove,idx)
        end
    end
    #println("before delete")
    #println("udegree1", udegree1, " ", length(udegree1))
    deleteat!(udegree1,remove)
    nbc1 :: Int64 = length(udegree1)
    deleteat!(udegree1count,remove)
    #println("after delete")
    #println("udegree1", udegree1, " ", length(udegree1))
    #println("udegree1count", udegree1count )
    #println("core2nodes", core2nodes)
    nbc2 :: Int64 = length(core2nodes)
    position = findin(edges,core2nodes)
    len = length(position)
    
    if core2nodes == unique(core2nodes)
        println("core2nodes unique!!")
    else
        println("WARNING!!")
        exit()
    end
    if udegree1 == unique(udegree1)
        println("udegree1 unique!!")
    else
        println("WARNING!!")
        exit()
    end
    if core1nodes == unique(core1nodes)
        println("core1nodes unique!!")
    else
        println("WARNING!!")
        exit()
    end
    #links = Array{Int64,1}(len)
    #linksneighbor = Array{Int64,1}(len)
    #### TODO: maybe this should be the length of unique core2nodes
    links = Array{Int64,1}()
    linksneighbor = Array{Int64,1}()
    #println("position length ", length(position), " should be bigger than the length core2nodes ", length(core2nodes))
    
    ### links, linksneighbor correspond to edges of core2node
    ### links and linksneighbor have the same size, but
    ### core2node is basically unique(links) and should have
    ### the same size as the first dimension of sizes.
    #println(position)
    for (idx, i) in enumerate(position)
        if isodd(i) == true
            if tms[edges[i]] > 1 && length(outneighbors(g, edges[i+1])) > 1
                push!(links, edges[i])
                push!(linksneighbor, edges[i + 1])
            end
        else
            if tms[edges[i]] > 1 && length(outneighbors(g, edges[i-1])) > 1
                push!(links, edges[i])
                push!(linksneighbor, edges[i - 1])
            end
        end
    end
    #println("links length ",length(links), " should be less than core2nodes length ", length(core2nodes))
    ### at this point we have:
    ### in core1nodes all nodes with tms 1.
    ### It means that we should remove some of them (those that are present in bridges).
    ### In core2nodes we have the opposite:
    ### we have all nodes that belong to a bigger component
    ### but we may need to add core2nodes that are part
    ### of smaller component (or even single-node component)
    ### but we don't know it yet.
    ###### TODO: check if order is kept!!
    ### The following line is an overestimation of the size of sizes
    # sizes = Array{Array{Int, 1}}(nbc2 + nbc1)
     sizes = Array{Array{Int, 1}}(nbc2)
    #sizes = Array{Array{Int, 1}}(length(links))
    for i in indices(sizes,1) sizes[i] = []; end

    bridges =  Array{Int64,1}()
#### TODO: maybe visited should be of size all edges of bridges type
#### to avoid duplicate calls to bfs_edge_subtree3!!
    visited = falses(length(links))
    boundary = Array{Int64,1}()

    println("Number of iterations: ", length(links))
    println("smaller number of iter ? : ", length(core2nodes))
    println("### setting up time:",time()- t, "(s)")
    ##### TODO: make sure that linksneighbor are unique
    ##### while links are not !!!!
    #println("links = ", links)
    #println("links length = ", length(links))
    #println("linksneighbor = ", linksneighbor)
    
    tr1 = time()
    for (i, u) in enumerate(links)
        v = linksneighbor[i]    
        ### this is directly a bridge
        if tms[v] > 1
            ## avoid double entries in bridges
            if (getindex(visited[i]) == false &&
                getindex(visited[findin(links, v)]) == false)
                push!(bridges, u)
                visited[i] = true
                push!(bridges, v)
                visited[findin(links, v)] = true
            end
            ### this is a degree one node
        ##### TODO: the following code should be deleted!!!    
        elseif length(outneighbors(g, v)) == 1
            println("WARNING: program should never reach these lines!! ")
            exit()
            unii = getindex(findin(core2nodes, u))
            len = length(sizes[unii])
            push!(boundary, u)
            if len == 0
                push!(sizes[unii], 1)
            else
                sizes[unii][1] += 1
            end
        else
            tree, distances, extra = bfs_edge_subtree3( g, u, v, tms)
            if isempty(extra) == true
                unii = getindex(findin(core2nodes, u))
                len = length(sizes[unii])
                cntr = zeros(Int64, maximum(distances))
                cntr = count3(distances,length(distances))
                #println(u, ",",v)
                #print(" $u with sizes[$unii]",sizes[unii] )
                #### code to have unique boundary nodes!
                if len == 0
                    #print(" ... will be pushed to boundary")
                    push!(boundary, u)
                else
                    #print(" (already in boundary)")
                end
                for (j,c) in enumerate(cntr)
                    if j <= len 
                        sizes[unii][j] += c
                    else
                        push!(sizes[unii], c)
                    end 
                end
                #println(" ... and after: ",sizes[unii] )
            else  ###
                #### means that we have intermediate nodes to be single node components!
                #### we make them all singular components
                #### TODO : I need a visited here!!!
                for i in 1:length(extra)
                    if extra[i]!=0
                        push!(bridges, i)
                        push!(bridges, extra[i])
                    end
                end
                push!(bridges, u)
                push!(bridges, tree[1])
                
            end
        end
    end
println("### exploration time:",time()- tr1, "(s)")

tr3 = time()
#println(core1nodes)
### make sure that bridges are not part of core1nodes
deleteat!(core1nodes,findin(core1nodes, bridges))

if boundary == unique(boundary)
    println("boundary unique!!")
else
    println("WARNING!!")
    exit()
end


additional = setdiff(udegree1, boundary)
#println("additional = ", additional)
#println("udegree1 = ", udegree1)
#println(length(boundary))
#println(length(additional))
#println("core2nodes = ", core2nodes)
new_sizes = Array{Array{Int64,1}}(length(boundary) + length(additional))
for i in indices(new_sizes,1) new_sizes[i] = []; end


j :: Int64 = 1
for i in 1:length(sizes)
    if !isempty(sizes[i])
        new_sizes[j] = sizes[i]
        #println("j =$j, i =$i ", core2nodes[i])
        u = core2nodes[i]
        if in(u, udegree1) == true 
            #println("new_sizes[$i] of $u should be changed: ", new_sizes[j][1])
            #println(findin(udegree1, u), " ", udegree1count[findin(udegree1, u)] )
            new_sizes[j][1] += getindex(udegree1count[findin(udegree1, u)])
        end
        j += 1
    end
end

## debugging
# println("new_sizes: ", new_sizes, " len = ", length(new_sizes), " j = $j")
# for i in 1:length(new_sizes)
#     if isempty(new_sizes[i]) == true
#         println("j must be: ", i)
#         break
#     end
# end

### updating additional nodes, that only have degree1 connections
for (i, u) in enumerate(udegree1)
    if in(u, additional) == true
        new_sizes[j] = [udegree1count[i]]
        j += 1
    end
end

boundary = [boundary ; additional]

#println("bridges = ",bridges)
#println("bdry = ", boundary)    
#println("sizes = ",new_sizes)
#println("core1nodes = ",core1nodes)



if length(new_sizes) != length(boundary)
    println("WARNING: Wrong size of sizes and core2nodes!")
    println("WARNING: Program will exit!")
    exit()
end
for u in core1nodes
    if in(u,boundary) == true || in(u,bridges) == true
        println("WARNING: core1nodes are not correct!")
        exit()
    end
end

println("### reshaping bridges time:",time()- tr3, "(s)")

B = Bridges(bridges, boundary, new_sizes)
println("### cleaning up time:",time()- t, "(s)")
return B, core1nodes
end


