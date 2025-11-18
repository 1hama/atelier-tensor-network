using ITensors, Plots, StatsBase, Luxor
include("../ctmrg/CTMRG_main.jl")

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

function ctmrg_snapshot(;W, boundary=nothing, maxsize, chi, width, height)
    #=
    W       : SYMMETRICAL ITensor representing the tensor network 
    maxsize : maximum number of iterations
    chi     : maximum bond dimension
    =#

    i, j, k, l = inds(W)
    W_bond_dim = dim(i)

    if isnothing(boundary)
        ones = ITensor(l)
        for val in 1:dim(l)
            ones[l=>val] = 1.0
        end
        boundary = ones
    end

    C, P, m, n, l = doCTMRG(;W=W, boundary=ones, maxsize=maxsize, chi=chi)

    # Snapshot in a row
    #=
                  C-m
      Ln[1] = L = P-s
                  C-n
    
                     P P P ... P-m
      L[width+1] = L W W W ... W-s
                     P P P ... P-n
                     #ofWs = width
    
                        m-C
      Rn[width+1] = R = s-P
                        n-C
    
              m-P P P ... P
      Rn[1] = s-W W W ... W R
              n-P P P ... P
              #ofWs = width
    =#

    s = Index(W_bond_dim, "s")
    L = replaceind(C * prime(C, m), m', n) # L has indices m, n, s
    L = prime(C,m) * replaceinds(P, [m, n, l]=>[m', n', s]) * prime(C, n)
    R = copy(L) # R has indices m, n, s

    # R[end] is equal to R so there has to be width + 1 elements in Rn
    Rn = [ITensor() for _ in 1:width+1]
    Ln = [ITensor() for _ in 1:width+1]
    Rn[end]   = copy(R)
    Ln[begin] = copy(L) # Initialize Ln from left to right
    
    # Store fixed Ws
    # (1,1) (1,2) ...
    # (2,1) (2,2) ...
    #  ...   ...  ...
    fixed_Ws       = Array{ITensor}(undef, height, width)
    fixed_Ws_ijkl  = Array{Tuple  }(undef, height, width)
    fixed_Ws_value = Array{Float64}(undef, height, width)

    # Initialize Rn from right to left
    for col in width:-1:1
        Rn[col] = (replaceinds(Rn[col + 1], [m, n, s]=>[m', n', s']) * prime(P, m)) *
                replaceinds(W, [i, j, k]=>[s, l', s']) * replaceinds(P, [n, l]=>[n', l'])
        Rn[col] /= norm(Rn[col]) # normalization for preventing overflow
    end
    @assert all(Rn .!= 0) "Rn still contains zero elements after initialization and update."

    function stochastic_index_sampling(environment)
        weights = []
        indices = Tuple{Int,Int,Int,Int}[]
        for (i_tmp, j_tmp, k_tmp, l_tmp) in Iterators.product(1:W_bond_dim, 1:W_bond_dim, 1:W_bond_dim, 1:W_bond_dim)
            weight = environment[i=>i_tmp, j=>j_tmp, k=>k_tmp, l=>l_tmp] * W[i=>i_tmp, j=>j_tmp, k=>k_tmp, l=>l_tmp]
            push!(weights, weight)
            push!(indices, (i_tmp, j_tmp, k_tmp, l_tmp))
        end

        sampled_idx = sample(Weights(weights ./ sum(weights)))
        sampled_ijkl = indices[sampled_idx]
        
        W_value = W[i=>sampled_ijkl[1], j=>sampled_ijkl[2], k=>sampled_ijkl[3], l=>sampled_ijkl[4]]
        fixed_W = W_value * onehot(i=>sampled_ijkl[1], j=>sampled_ijkl[2], k=>sampled_ijkl[3], l=>sampled_ijkl[4])
        return fixed_W, sampled_ijkl, W_value
    end

    # Fix Ws of the first row from left to right stochastically
    for col in 1:width
        # environment: i,j,k,l as indices
        environment = replaceind(Ln[col], s=>i) * (replaceinds(Rn[col + 1], [m, n, s]=>[m', n', k]) *
                prime(P, m)) * replaceinds(P, [n, l]=>[n', j])

        fixed_Ws[1,col], fixed_Ws_ijkl[1,col], fixed_Ws_value[1,col] = stochastic_index_sampling(environment)

        Ln[col + 1] = replaceinds(replaceinds(Ln[col], [m, n, s]=>[m', n', s']) * replaceinds(P, [m, l]=>[m', j]) *
                replaceinds(fixed_Ws[1,col], [i, k]=>[s', s]) * replaceinds(P, [n, l]=>[n', l]), [m, n]=>[n, m])
        Ln[col + 1] /= norm(Ln[col + 1]) # normalization for preventing overflow
    end
    @assert all(Ln .!= 0) "Ln still contains zero elements after initialization and update."

    top_left = copy(C)
    top_right = copy(C)
    for row in 2:height
        top_left  = replaceind(top_left  * prime(P, m) * onehot(l=>fixed_Ws_ijkl[row-1, 1][1]    ), m'=>n)
        top_right = replaceind(top_right * prime(P, m) * onehot(l=>fixed_Ws_ijkl[row-1, width][3]), m'=>n)
        L = replaceind(top_left  * prime(P, m) * prime(C, m), l=>s)
        R = replaceind(top_right * prime(P, m) * prime(C, m), l=>s)
        
        Rn = [ITensor() for _ in 1:width+1]
        Ln = [ITensor() for _ in 1:width+1]
        Rn[end]   = R
        Ln[begin] = L
        
        # initialize Rn
        for col in width:-1:1
            Rn[col] = replaceinds(replaceinds(Rn[col + 1], [s, n]=>[k, n']) * P *
                    (fixed_Ws_value[row-1, col] * onehot(l=>fixed_Ws_ijkl[row-1, col][2], j=>fixed_Ws_ijkl[row-1, col][4])) *
                    W * prime(P, n), [m, n, i]=>[n, m, s])
            Rn[col] /= norm(Rn[col]) # normalization for preventing overflow
        end

        # fix spins stochastically
        for col in 1:width
            # environment: i,j,k,l as indices
            environment = replaceind(Ln[col], s=>i) *
                    (prime(P, n) * fixed_Ws_value[row-1, col] *
                    onehot(l=>fixed_Ws_ijkl[row-1, col][2], j=>fixed_Ws_ijkl[row-1, col][4])) *
                    P * replaceinds(Rn[col + 1], [m, n, s]=>[n', m, k])

            fixed_Ws[row,col], fixed_Ws_ijkl[row,col], fixed_Ws_value[row,col] = stochastic_index_sampling(environment)

            Ln[col + 1] = replaceinds(replaceinds(Ln[col], [s, n]=>[i, n']) * P *
                    (fixed_Ws_value[row-1, col] * onehot(l=>fixed_Ws_ijkl[row-1, col][2], j=>fixed_Ws_ijkl[row-1, col][4])) *
                    fixed_Ws[row,col] * prime(P, n), [m, n, k]=>[n, m, s])
            Ln[col + 1] /= norm(Ln[col + 1]) # normalization for preventing overflow
        end
    end

    return fixed_Ws_ijkl
end

function get_snapshot_output_path(base_dir, filename)
    dir = joinpath(dirname(base_dir), "snapshot_output")
    mkpath(dir)
    return joinpath(dir, filename)
end