using Laplacians

function LinvdiagSS(a::SparseMatrixCSC{Float64}; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

    f = approxCholLap(a,tol=1e-5);
    
    n = size(a,1)
    k = round(Int, JLfac*log(n)) # number of dims for JL
    
    U = wtedEdgeVertexMat(a)
    m = size(U,1)
    er = zeros(n)

    # v_i
    for i = 1:k # q 
        # random gaussian projections Q
        r = randn(m) 
        # compute (QW^(1/2)B)
        ur = U'*r 
        v = zeros(n)
        # solve(L, ((QW^(1/2)B)^T))^T
        # v here is Z_{i,:}
        v = f(ur[:])
        # er(u) = Sigma_{v\in V\u}(Z_{u,v})^2
        # decause I'll do the above k times
        # given the k random projections,
        # I have to divide by k
        er.+= v.^2/k
    end

  return er;

end
