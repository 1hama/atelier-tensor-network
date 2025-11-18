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

    
    # Spin-spin correlation
    SS = ITensor(i, j, k, l)
    SS[i=>1, j=>1, k=>1, l=>1] = W[i=>1, j=>1, k=>1, l=>1] * 1
    SS[i=>2, j=>2, k=>2, l=>2] = W[i=>2, j=>2, k=>2, l=>2] * 1
    SS[i=>1, j=>2, k=>1, l=>2] = W[i=>1, j=>2, k=>1, l=>2] * (-1)
    SS[i=>2, j=>1, k=>2, l=>1] = W[i=>2, j=>1, k=>2, l=>1] * (-1)
    spin_spin_col = only(environment_tensor * SS / (environment_tensor * W))
    return spin_spin_col
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

# exact expression for nearest-neighbor spin correlation
function nn_corr(T)
    β = 1 / T
    k = 1 / sinh(2β)^2
    m = (2*sqrt(k)/(1+k))^2  # elliptic modulus squared
    Kval = ellipk(m)         # complete elliptic integral of the first kind
    return 0.5 * coth(2β) * (1 + (2/π) * (2*tanh(2β)^2 - 1) * Kval)
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

function nn_corr_deriv(x, dx=1e-5)
    (nn_corr(x+dx) - nn_corr(x-dx)) / (2*dx)
end

function main()
    # Define constants
    J = 1.0
    chi = 20
    maxsize = 100

    xlims = (1, 5)
    Tc = 2 / log(1+sqrt(2))
    # hs = [0, 1, 2, 3]
    hs = range(0, 0.4, length = 5)
    markers = [:circle, :cross, :utriangle, :square, :pentagon]
    colors = [:black, :black, :black, :black, :black]
    sizes = [2.5,2.5,2.5,2.5,2.5,]

    set_Plots_default()
    plt = plot(;
        xlims = xlims,
        ylims = (0.2, 1.0),
        xlabel = "\$T\$",
        xguidefont = font(22),
        xtickfont = font(12),
        ytickfont = font(12),
        # title = "J = $J, chi = $chi, maxsize = $maxsize",
        legend = :topright,
        legendfont = 10
    )

    # Draw theoretical curve for h = 0 only
    plot!(x -> nn_corr(x), xlims[1], xlims[2]; color = :gray, lw = 0.8, seriestype = :line, marker = :none, label = "Theory ( \$h=0\$ )")

    # Plot computed values for each h
    for (idx, h) in enumerate(hs)
        xs = sample_by_density(xlims[1], xlims[2], 40, x -> atan(2 * nn_corr_deriv(x)))
        # xs = range(xlims[1], xlims[2], length = 40)
        results = Vector{Float64}(undef, length(xs))
        @showprogress 1 "Running CTMRG for h=$h..." for (i, x) in enumerate(xs)
            results[i] = ctmrg(
                W = create_W(; beta = 1 / x, J = J, h = h),
                maxsize = maxsize,
                chi = chi
            )
        end
        plot!(
            xs, results;
            color = colors[idx],
            marker = markers[idx],
            markersize = sizes[idx],
            label = "\$ h = $h \$",
            seriestype = :scatter
        )
    end

    # Draw vertical line at Tc
    vline!([Tc], color = :black, linestyle = :dash, lw = 1, label="")

    annotate!(Tc, 0.14, Plots.text("\$T_{\\mathrm{c}}\$", 18))

    filename = "ising_multi_h_chi-$(chi)_N-$(maxsize)"
    output_path = get_graph_output_path(@__DIR__, filename * ".pdf")
    savefig(output_path)
    println("SVG image saved at $(output_path)")
    display(plt)
end

main()