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
                        #push!(core2nodes, w)
                        ########
                        #edge = Edge(v,w)
                        #push!(bridges, edge)
                    elseif length(outneighbors(g, w)) == 1
                        push!(core1nodes, w)
                        push!(core1neighbor, v)
                        #push!(core2nodes, v)
                        #########
                        #edge = Edge(v,w)
                        #push!(bridges, edge)
                    ### I do not care if edges are ordered    
                    # elseif v < w
                    #     edge = Edge(v, w);
                    #     push!(core3nodes,v);
                    #     push!(core3nodes,w);
                    #     push!(bridges, edge)
                        
                    else
                        ### I do not care about order!!
                        # edge = Edge(w, v)
                        # push!(core3nodes,w);
                        # push!(core3nodes,v);
                        # push!(bridges, edge)
                        push!(edges,w)
                        push!(edges,v)
                    end
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
    ## core2pairs is replaced by sizes
    # core2pairs = zeros(Int64, length(core2nodes))
    # for (idx, u) in enumerate(core2nodes)
    #     c = count(x->x==u,core1neighbor);
    #     core2pairs[idx] = c;
    # end
    # println(core2pairs)
    t = time()
    ###### TODO: sort(unique()) or unique(sort()) ?
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
    B = Bridges(edges, Set(core2nodes), sizes)
    ### remove the creation of core3nodes!!! and call following constructor
    #B = Bridges(edges, core2nodes, core2pairs)
    
    return B, core1nodes
end

