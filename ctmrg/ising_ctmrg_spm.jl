using ITensors, Plots, SpecialFunctions, ProgressMeter
include("CTMRG_main.jl")

#=
 j
iWk
 l

boundary-l
=#

function ctmrg(;W, maxsize=30, chi=10)
    #=
    W       : SYMMETRIC ITensor representing the tensor network 
    maxsize : maximum number of iterations
    chi     : maximum bond dimension
    =#
    i, j, k, l = inds(W)
    W_bond_dim = dim(i)

    V = ITensor(l)
    V[l=>1] = sqrt(W[i=>1, j=>1, k=>1, l=>2])
    V[l=>2] = sqrt(W[i=>1, j=>2, k=>2, l=>2])

    C, P, m, n, l = doCTMRG(;W=W, boundary=V, maxsize=maxsize, chi=chi)

    # C index m, n
    # P index m, n, l

    bottom_half = C * prime(P, n) * prime(C, n)

    environment_tensor = replaceinds(bottom_half, [n, m]=>[n', m']) * replaceinds(P, [m, l]=>[m', i]) * 
            replaceinds(P, [n, l]=>[n', k]) * replaceind(bottom_half, l=>j)

    # Spin average
    S = ITensor(i, j, k, l)
    S[i=>1, j=>1, k=>1, l=>1] = W[i=>1, j=>1, k=>1, l=>1] * 4
    add_symmetrically(S, (i=>2, j=>1, k=>1, l=>1), W[i=>2, j=>1, k=>1, l=>1] * 3)
    add_symmetrically(S, (i=>2, j=>1, k=>2, l=>1), W[i=>2, j=>1, k=>2, l=>1] * 2)
    add_symmetrically(S, (i=>2, j=>2, k=>1, l=>1), W[i=>2, j=>2, k=>1, l=>1] * 2)
    add_symmetrically(S, (i=>1, j=>2, k=>2, l=>2), W[i=>1, j=>2, k=>2, l=>2] * 1)
    S[i=>2, j=>2, k=>2, l=>2] = W[i=>2, j=>2, k=>2, l=>2] * 0
    spin_average = only(environment_tensor * S / (environment_tensor * W)) / 4
    return spin_average
end

function create_W(;beta, J, h)
    i = Index(2, "i")
    j = Index(2, "j")
    k = Index(2, "k")
    l = Index(2, "l")
    W = ITensor(Float64, i, j, k, l)

    # 1 for ↑ and 2 for ↓
    W[i=>1, j=>1, k=>1, l=>1] = exp(beta*(4*J + 2*h))
    add_symmetrically(W, (i=>1, j=>1, k=>1, l=>2), exp(beta*h))
    add_symmetrically(W, (i=>1, j=>1, k=>2, l=>2), 1)
    add_symmetrically(W, (i=>1, j=>2, k=>1, l=>2), exp(-4*beta*J))
    add_symmetrically(W, (i=>1, j=>2, k=>2, l=>2), exp(-beta*h))
    W[i=>2, j=>2, k=>2, l=>2] = exp(beta*(4*J - 2*h))

    return W
end

function sample_by_density(xmin, xmax, npoints, rho)
    """
    sample values with density rho
    """
    function interp1(cdf, xgrid, p)
        idx = findfirst(c -> c >= p, cdf)
        if idx == 1 || idx === nothing
            return xgrid[1]
        elseif idx == length(xgrid)
            return xgrid[end]
        else
            # Linear interpolation between xgrid[idx-1] and xgrid[idx]
            c1, c2 = cdf[idx-1], cdf[idx]
            x1, x2 = xgrid[idx-1], xgrid[idx]
            return x1 + (p - c1) * (x2 - x1) / (c2 - c1)
        end
    end
    # Create a fine grid for integration
    xgrid = range(xmin, xmax, length=10000)
    dx = (xmax - xmin) / (length(xgrid) - 1)
    # Compute cumulative density (CDF)
    cdf = cumsum([rho(x)*dx for x in xgrid])
    cdf ./= cdf[end]  # Normalize to [0,1]
    # Invert CDF to get xs
    xs = [interp1(cdf, xgrid, p) for p in range(0, 1, length=npoints)]
    return xs
end

function main()
    # Define constants
    J = 1.0
    # k_B = 1.0 : Boltzmann constant
    chi = 20
    maxsize = 100

    xmax = 4
    hmin = -0.005
    hmax = 0.005

    # critical temperature
    Tc = 2 / log(1+sqrt(2))

    default(legend = false)
    xs = sample_by_density(1, xmax, 20, x -> exp(-0.5 * ((x - Tc)/0.4)^2))
    hs = range(hmin, hmax, length = 21)

    # Precompute zvals
    zvals = Array{Float64}(undef, length(hs), length(xs))

    println("Starting simulation...")
    start_time = time()
    @showprogress 1 "Running CTMRG..." for (j, hval) in enumerate(hs)
        for (i, xval) in enumerate(xs)
            zvals[j, i] = ctmrg(
                W = create_W(; beta = 1 / xval, J = J, h = hval),
                maxsize = maxsize,
                chi = chi
            )
        end
    end
    elapsed = time() - start_time
    println("Execution time: $(round(elapsed, digits=2)) seconds")

    p = plot(
        xs, hs, zvals,
        st = :wireframe,    # wireframe grid only
        c = :black,         # black lines for print
        linewidth = 0.5,
        grid = true,
        gridalpha = 1,
        gridlinewidth = 0.5,
        xflip = true,
        yflip = false,
        xlims = (1, xmax),     # set display range for x
        ylims = (hmin, hmax),    # set display range for h
        zlims = (0, 1)
    )

    plot!(p, camera = (20, 40))

    filename = "ising_spm_chi-$(chi)_N-$(maxsize)"
    output_path = get_graph_output_path(@__DIR__, filename * ".pdf")
    savefig(output_path)
    println("SVG image saved at $(output_path)")
    display(p)
end

main()