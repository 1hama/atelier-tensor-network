using ITensors, Plots
include("CTMRG_main.jl")

#=
 j
iWk
 l
=#

#=
 j      l     m
iWk -> mPn , nC
 l

ones-l
=#

function ctmrg(;W, maxsize=30, chi=10)
    #=
    W       : SYMMETRIC ITensor representing the tensor network 
    maxsize : maximum number of iterations
    chi     : maximum bond dimension
    =#
    i, j, k, l = inds(W)
    W_bond_dim = dim(i)

    ones = ITensor(l)
    for val in 1:dim(l)
        ones[l=>val] = 1.0
    end

    C, P, m, n, l = doCTMRG(;W=W, boundary=ones, maxsize=maxsize, chi=chi)

    # C index m, n
    # P index m, n, l

    bottom_half = C * prime(P, n) * prime(C, n)

    # ======METHOD 1======
    # spin = ITensor(l, l')
    # spin[l=>1, l'=>1] = -1.0
    # spin[l=>2, l'=>2] = 1.0
    # spin[l=>1, l'=>2] = 0.0
    # spin[l=>2, l'=>1] = 0.0

    # spin_average = only((prime(bottom_half,l) * spin * bottom_half) / (bottom_half * bottom_half))

    # return spin_average

    # ======METHOD 2======

    environment_tensor = replaceinds(bottom_half, [n, m]=>[n', m']) * replaceinds(P, [m, l]=>[m', i]) * 
            replaceinds(P, [n, l]=>[n', k]) * replaceind(bottom_half, l=>j)

    SS = ITensor(i, j, k, l)
    SS[i=>1, j=>1, k=>1, l=>1] = W[i=>1, j=>1, k=>1, l=>1] * 1
    SS[i=>2, j=>2, k=>2, l=>2] = W[i=>2, j=>2, k=>2, l=>2] * 1
    SS[i=>1, j=>2, k=>1, l=>2] = W[i=>1, j=>2, k=>1, l=>2] * (-1)
    SS[i=>2, j=>1, k=>2, l=>1] = W[i=>2, j=>1, k=>2, l=>1] * (-1)

    spin_spin_col = only(environment_tensor * SS / (environment_tensor * W))

    return spin_spin_col
end

function create_W(x, a)
    i = Index(2, "i")
    j = Index(2, "j")
    k = Index(2, "k")
    l = Index(2, "l")
    W = ITensor(i, j, k, l)

    # 1 for ↑↓ and ↓↑, 2 for ↑↑ and ↓↓
    W[i=>1, j=>1, k=>1, l=>1] = exp(4*a/x)
    W[i=>2, j=>2, k=>2, l=>2] = exp(4*a/x)
    add_symmetrically(W, (i=>1, j=>2, k=>1, l=>2), exp(-4*a/x))
    add_symmetrically(W, (i=>1, j=>1, k=>1, l=>2), 1)
    add_symmetrically(W, (i=>1, j=>1, k=>2, l=>2), 1)
    add_symmetrically(W, (i=>1, j=>2, k=>2, l=>2), 1)

    return W
end

function plot_function(f; xs=nothing, xlims=(1, 5), ylims=(0.2, 1.0), npoints=100, title="Function Plot", vline_x=nothing, filename)
    if xs === nothing
        xs = range(xlims[1], xlims[2], length=npoints)
    end
    ys = [f(x) for x in xs]

    set_Plots_default()
    plt = plot(xs, ys;
        xlims=xlims,
        ylims=ylims,
        xlabel="\$T\$",
    )
    if vline_x !== nothing
        vline!([vline_x], color=:black, linestyle=:dash)
        annotate!(vline_x, 0.15, Plots.text("\$T=$vline_x\$", 10))
    end
    output_path = get_graph_output_path(@__DIR__, filename*".pdf")
    savefig(output_path)
    println("SVG image saved at $(output_path)")
    display(plt)
end

function main()
    # Define constants
    J = 1.0 
    k_B = 1.0 # Boltzmann constant
    a = J / k_B
    chi = 40
    maxsize = 20

    # Example: denser sampling near x=2.27
    xs1 = range(1, 2.0, length=20)
    xs2 = range(2.0, 2.5, length=60)
    xs3 = range(2.5, 5, length=50)
    xs = vcat(collect(xs1), collect(xs2), collect(xs3))

    println("Starting simulation...")
    start_time = time()
    filename = "\$ ising_a-$(a)_chi-$(chi)_N-$(maxsize)_\$"
    title = "a = $a, chi = $chi, maxsize = $maxsize"
    plot_function(x -> ctmrg(; W=create_W(x, a), maxsize=maxsize, chi=chi), 
        xs=xs, title=title, vline_x=2.27, filename=filename)
    elapsed = time() - start_time
    println("Execution time: $(round(elapsed, digits=2)) seconds")
end

main()