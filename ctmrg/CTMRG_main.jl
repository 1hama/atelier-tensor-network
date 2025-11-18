using ITensors

function doCTMRG(;W, boundary, maxsize=30, chi=10)
    """
    W       : SYMMETRIC ITensor representing the tensor network 
    V       : boundary vector (ITensor)
    maxsize : maximum number of iterations
    chi     : maximum bond dimension
    """
    if length(inds(W)) != 4
        error("W must have exactly 4 indices")
    end
    if !isSymmetric(W)
        error("W must be symmetric")
    end
    if length(inds(boundary)) != 1
        error("Boundary must have only one index")
    end
    if !(inds(boundary)[1] == inds(W)[4])
        error("Boundary index must be the fourth indices of W")
    end

    #=
     j      l     m
    iWk -> mPn , nC
     l

    boundary-l
    =#

    i, j, k, l = inds(W)
    W_bond_dim = dim(i)
    # initiate P, C
    m = Index(W_bond_dim, "m")
    n = Index(W_bond_dim, "n")
    P = replaceinds(W * boundary, [i, j, k], [m, l, n])
    C = replaceinds(W * boundary * (boundary * delta(l, k)), [i, j], [n, m])

    for size in 1:maxsize-1
        # update C and P
        Pnext = P * W
        C = Pnext * prime(C, n) * replaceinds(P, [n, l]=>[n', i])
        A, S, A_daggar = svd(C, j, m; maxdim=chi)
        link1, link2 = inds(S)
        m = Index(dim(link1), "m")
        n = Index(dim(link1), "n")
        P = replaceinds(Pnext * replaceind(A, j=>i) * A_daggar, [link1, link2, j]=>[m, n, l])
        C = replaceinds(S, [link1, link2], [m, n])

        # normalize C and P
        C /= norm(C)
        P /= norm(P)
        # @show size, norm(C), norm(P)
    end

    # C index: m, n
    # P index: m, n, l

    return C, P, m, n, l
end

function isSymmetric(W)
    # Check if W is symmetric under rotation and reflection
    # W should be a 4-index ITensor with indices i, j, k, l
    i, j, k, l = inds(W)
    if length(inds(W)) != 4
        return false
    end

    # Check if W is symmetric under rotation (cycle indices)
    for a in 1:dim(i), b in 1:dim(j), c in 1:dim(k), d in 1:dim(l)
        v = (a, b, c, d)
        if !(W[i=>v[1], j=>v[2], k=>v[3], l=>v[4]] ==
            W[j=>v[1], k=>v[2], l=>v[3], i=>v[4]] ==
            W[k=>v[1], l=>v[2], i=>v[3], j=>v[4]] ==
            W[l=>v[1], i=>v[2], j=>v[3], k=>v[4]])
            println("Rotation symmetry failed for indices: ", v)
            return false
        end
    end

    # Check if W is symmetric under reflection (reverse all indices)
    for a in 1:dim(i), b in 1:dim(j), c in 1:dim(k), d in 1:dim(l)
        if W[i=>a, j=>b, k=>c, l=>d] != W[i=>d, j=>c, k=>b, l=>a]
            println("Reflection symmetry failed for indices: ", (a, b, c, d))
            return false
        end
    end
    return true
end

function add_symmetrically(W::ITensor, ijklvals::NTuple{4, Pair{<:Index, Int}}, val::Union{Float64, Int})
    # modify W (which is mutable) to be symmetrical under rotation and reflection

    i, j, k, l = [p.first for p in ijklvals]
    v1, v2, v3, v4 = [p.second for p in ijklvals]
    
    W[i=>v1, j=>v2, k=>v3, l=>v4] = val
    W[i=>v4, j=>v1, k=>v2, l=>v3] = val
    W[i=>v3, j=>v4, k=>v1, l=>v2] = val
    W[i=>v2, j=>v3, k=>v4, l=>v1] = val

    W[i=>v1, j=>v4, k=>v3, l=>v2] = val
    W[i=>v2, j=>v1, k=>v4, l=>v3] = val
    W[i=>v3, j=>v2, k=>v1, l=>v4] = val
    W[i=>v4, j=>v3, k=>v2, l=>v1] = val

    return W
end

function get_graph_output_path(base_dir, filename)
    dir = joinpath(dirname(base_dir), "graph_output")
    mkpath(dir)
    return joinpath(dir, filename)
end

function set_Plots_default()
    Plots.default(
        seriestype=:scatter,
        marker=:circle,
        color=:black,
        markersize=1.5,
        xguidefontsize=15,
        # ylabel="",
        # title=title,
        legend=false,
        xminorticks=false,
        grid=false,
        framestyle=:box,
    )
end