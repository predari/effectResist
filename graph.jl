struct Graph
    n :: Int # |V|
    m :: Int # |E|
    nbr :: Array{Array{Int, 1}, 1} # neighbros of each vertex
end

# Graph -> Graph
function max_cc(G :: Graph)
    vis = fill(false, G.n)

    dfs(u :: Int) = begin
        vis[u] = true
        ret = [u]
        for v in G.nbr[u]
            if ! vis[v]
                append!(ret, dfs(v))
            end
        end
        return ret
    end

    max_cc = Vector{Int}(0)
    for u in 1 : G.n
        if ! vis[u]
            cc = dfs(u)
            if length(cc) > length(max_cc)
                max_cc = cc
            end
        end
    end

    id = Dict{Int,Int}( tuple.(max_cc, indices(max_cc,1)) )

    n = length(id)

    nbr = Array{Array{Int, 1}}(n)
    for u in filter( x -> haskey(id, x), indices(G.nbr, 1) )
        nbr[ id[u] ] = getindex.(id, filter( x -> haskey(id, x), G.nbr[u] ) )
    end

    m = div( sum( length.(nbr) ), 2 )

    return Graph(n, m, nbr)
end

# Graph -> Graph
function bfs_max_cc(G :: Graph)
    vis = fill(false, G.n)

    bfs(u :: Int) = begin
        vis[u] = true
        ret = [u]
        hinge = [u]
        while !isempty(hinge)
            v = pop!(hinge)
            for t in G.nbr[v]
                if !vis[t]
                    push!(ret, t)
                    push!(hinge,t)
                    vis[t] = true
                    #println(t)
                end
            end
        end
        return ret
    end
    max_cc = Vector{Int}(0)
    for u in 1 : G.n
        if ! vis[u]
            cc = bfs(u)
            if length(cc) > length(max_cc)
                max_cc = cc
            end
        end
    end

    id = Dict{Int,Int}( tuple.(max_cc, indices(max_cc,1)) )

    n = length(id)

    nbr = Array{Array{Int, 1}}(n)
    for u in filter( x -> haskey(id, x), indices(G.nbr, 1) )
        nbr[ id[u] ] = getindex.(id, filter( x -> haskey(id, x), G.nbr[u] ) )
    end

    m = div( sum( length.(nbr) ), 2 )

    return Graph(n, m, nbr)
end

# String -> Graph
function read_graph(str :: String)
    ints = parse.(Int, split(str))

    n = 0
    d = Dict{Int,Int}()
    edges = Set{Tuple{Int,Int}}()
    getid(x :: Int) = haskey(d, x) ? d[x] : d[x] = n += 1

    for i in 1 : 2 : length(ints)
        u = getid(ints[i])
        v = getid(ints[i + 1])
        if u == v continue end
        if u > v  u,v = v,u end
        push!(edges, (u,v))
    end
    m = length(edges)

    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    for (u,v) in edges
        push!(nbr[u], v)
        push!(nbr[v], u)
    end

    return Graph(n,m,nbr)
end

# String -> Graph
function read_file(filename :: String)
    return open(filename) do f
        #read_graph( readstring(f) )
        bfs_max_cc( read_graph( read(f, String) ) )
    end
end

function sparsemat2Graph(A :: SparseMatrixCSC{Float64})
    nbr = Array{Array{Int, 1}}(A.n)
    (I, J)=findnz(A)
    print("A.m=",A.m, " A.n=", A.n, " nnz(A)=",nnz(A), "\n")
    for u in 1:A.m
        st = 0
        fs = 0
        for idx in find(x->x==u, J)
            if st == 0
                st = idx
                fs = st
            else
                fs = fs + 1
            end
        end
        nbr[u] = I[st:fs]
    end
    return Graph(A.n, nnz(A), nbr)
end


function Components3Graph(n:: Int, m:: Int)
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 2)
    push!(nbr[2], 1)
    push!(nbr[1], 4)
    push!(nbr[4], 1)
    push!(nbr[2], 3)
    push!(nbr[3], 2)
    push!(nbr[2], 4)
    push!(nbr[4], 2)
    push!(nbr[3], 4)
    push!(nbr[4], 3)
    push!(nbr[4], 5)
    push!(nbr[5], 4)
    push!(nbr[5], 6)
    push!(nbr[6], 5)
    push!(nbr[5], 8)
    push!(nbr[8], 5)
    push!(nbr[6], 8)
    push!(nbr[8], 6)
    push!(nbr[6], 7)
    push!(nbr[7], 6)    
    return Graph(n,m,nbr)
end

function ComponentExtnodes3Graph(n:: Int, m:: Int)
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 2)
    push!(nbr[1], 7)
    push!(nbr[1], 3)
    push!(nbr[2], 1)
    push!(nbr[2], 4)
    push!(nbr[2], 3)
    push!(nbr[3], 1)
    push!(nbr[3], 2)
    push!(nbr[3], 4)
    push!(nbr[3], 5)
    push!(nbr[4], 3)
    push!(nbr[4], 2)
    push!(nbr[4], 6)
    push!(nbr[5], 3)
    push!(nbr[6], 4)
    push!(nbr[7], 1)
    return Graph(n,m,nbr)
end

function Component2sExtnodes3Graph(n:: Int, m:: Int)
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 2)
    push!(nbr[1], 7)
    push!(nbr[1], 3)
    push!(nbr[2], 1)
    push!(nbr[2], 4)
    push!(nbr[2], 3)
    push!(nbr[2], 8)
    push!(nbr[3], 1)
    push!(nbr[3], 2)
    push!(nbr[3], 4)
    push!(nbr[3], 5)
    push!(nbr[4], 3)
    push!(nbr[4], 2)
    push!(nbr[4], 6)
    push!(nbr[5], 3)
    push!(nbr[6], 4)
    push!(nbr[7], 1)
    push!(nbr[8], 2)
    push!(nbr[8], 9)
    push!(nbr[8], 10)
    push!(nbr[9], 8)
    push!(nbr[9], 10)
    push!(nbr[10], 8)
    push!(nbr[10], 9)
    return Graph(n,m,nbr)
end

function LineGraph(n:: Int)
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    if n >= 2
        push!(nbr[1], 2)
    end
    for i in 2:n-1
        push!(nbr[i], i+1)
        push!(nbr[i], i-1)
    end
    push!(nbr[n], n-1)
    m = 2*n
    return Graph(n,m,nbr)
end


function AlmostLineGraph(n:: Int, m:: Int)
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 2)
    push!(nbr[2], 1)
    push!(nbr[2], 3)
    push!(nbr[2], 4)
    push!(nbr[3], 2)
    push!(nbr[4], 3)
    # line graph plus (2,3) edge
    push!(nbr[3], 4)
    push!(nbr[4], 2)
    return Graph(n,m,nbr)
end


function TriangleGraph(n:: Int, m:: Int)
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 2)
    push!(nbr[1], 3)
    push!(nbr[2], 3)
    push!(nbr[2], 1)
    push!(nbr[3], 1)
    push!(nbr[3], 2)
    return Graph(n,m,nbr)
end


function TestGraph(n:: Int, m:: Int)
    ## n = 15 m = 20*2
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 14)
    push!(nbr[1], 15)
    push!(nbr[1], 2)
    push!(nbr[1], 3)
    push!(nbr[1], 4)
    push!(nbr[2], 1)
    push!(nbr[2], 5)
    push!(nbr[2], 10)
    push!(nbr[2], 9)
    push!(nbr[3], 1)
    push!(nbr[3], 4)
    push!(nbr[3], 5)
    push!(nbr[4], 1)
    push!(nbr[4], 3)
    push!(nbr[4], 6)
    push!(nbr[4], 7)
    push!(nbr[5], 3)
    push!(nbr[5], 2)
    push!(nbr[5], 6)
    push!(nbr[5], 9)
    push!(nbr[6], 4)
    push!(nbr[6], 5)
    push!(nbr[6], 7)
    push!(nbr[7], 6)
    push!(nbr[7], 4)
    push!(nbr[7], 8)
    push!(nbr[8], 7)
    push!(nbr[9], 2)
    push!(nbr[9], 5)
    push!(nbr[10], 2)
    push!(nbr[10], 12)
    push!(nbr[10], 11)    
    push!(nbr[11], 12)
    push!(nbr[11], 10)
    push!(nbr[12], 11)
    push!(nbr[12], 10)
    push!(nbr[12], 13)
    push!(nbr[13], 12)
    push!(nbr[14], 1)
    push!(nbr[15], 1)
    
    return Graph(n,m,nbr)
end


function subTestGraph(n:: Int, m:: Int)
    ## n = 11 m = 15*2
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 10)
    push!(nbr[1], 11)
    push!(nbr[1], 2)
    push!(nbr[1], 3)
    push!(nbr[1], 4)
    push!(nbr[2], 1)
    push!(nbr[2], 5)
    push!(nbr[2], 9)
    push!(nbr[3], 1)
    push!(nbr[3], 4)
    push!(nbr[3], 5)
    push!(nbr[4], 1)
    push!(nbr[4], 3)
    push!(nbr[4], 6)
    push!(nbr[4], 7)
    push!(nbr[5], 3)
    push!(nbr[5], 2)
    push!(nbr[5], 6)
    push!(nbr[5], 9)
    push!(nbr[6], 4)
    push!(nbr[6], 5)
    push!(nbr[6], 7)
    push!(nbr[7], 6)
    push!(nbr[7], 4)
    push!(nbr[7], 8)
    push!(nbr[8], 7)
    push!(nbr[9], 2)
    push!(nbr[9], 5)
    push!(nbr[10], 1)
    push!(nbr[11], 1)
    
    return Graph(n,m,nbr)
end


function subTestGraph2(n:: Int, m:: Int)
    ## n = 8 m = 12*2
    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    push!(nbr[1], 2)
    push!(nbr[1], 3)
    push!(nbr[1], 4)
    push!(nbr[2], 1)
    push!(nbr[2], 5)
    push!(nbr[2], 8)
    push!(nbr[3], 1)
    push!(nbr[3], 4)
    push!(nbr[3], 5)
    push!(nbr[4], 1)
    push!(nbr[4], 3)
    push!(nbr[4], 6)
    push!(nbr[4], 7)
    push!(nbr[5], 3)
    push!(nbr[5], 2)
    push!(nbr[5], 6)
    push!(nbr[5], 8)
    push!(nbr[6], 4)
    push!(nbr[6], 5)
    push!(nbr[6], 7)
    push!(nbr[7], 6)
    push!(nbr[7], 4)
    push!(nbr[8], 2)
    push!(nbr[8], 5)
    
    return Graph(n,m,nbr)
end
