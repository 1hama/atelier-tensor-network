using ITensors, Plots, StatsBase, Luxor
include("snapshot.jl")

function draw_something(fixed_ijkl; cellsize=20, filename="something")
    # This function will draw (width-2) by (height-2) diagram

    height, width = size(fixed_ijkl)

    output_path = get_snapshot_output_path(@__DIR__, filename*".svg")

    Drawing((width-2) * cellsize, (height-2) * cellsize, output_path)
    origin(Point(-cellsize, -cellsize)) # origin at top left
    background("white")

    sethue("black")

    legwidth = cellsize / 4

    roundfactor = cellsize / 10

    function drawleg(x, y, n, action)
        if n == 1
            box(Point(x+cellsize/4+legwidth/4-roundfactor/2, y+cellsize/2), cellsize/2+legwidth/2+roundfactor, legwidth, roundfactor, action = action)
        elseif n == 2
            box(Point(x+cellsize/2, y+cellsize/4+legwidth/4-roundfactor/2), legwidth, cellsize/2+legwidth/2+roundfactor, roundfactor, action = action)
        elseif n == 3
            box(Point(x+cellsize/4*3-legwidth/4+roundfactor/2, y+cellsize/2), cellsize/2+legwidth/2+roundfactor, legwidth, roundfactor, action = action)
        elseif n == 4
            box(Point(x+cellsize/2, y+cellsize/4*3-legwidth/4+roundfactor/2), legwidth, cellsize/2+legwidth/2+roundfactor, roundfactor, action = action)
        end
    end

    for row in 2:height-1
        for col in 2:width-1
            x = (col - 1) * cellsize
            y = (row - 1) * cellsize
            ijkl = fixed_ijkl[row, col]

            for leg in 1:4
                if ijkl[leg] == 2
                    sethue("gray")
                    legwidth = cellsize/6
                    drawleg(x, y, leg, :fill)
                end
            end
        end
    end
    for row in 2:height-1
        for col in 2:width-1
            x = (col - 1) * cellsize
            y = (row - 1) * cellsize
            ijkl = fixed_ijkl[row, col]

            for leg in 1:4
                if ijkl[leg] == 3
                    sethue("black")
                    legwidth = cellsize/4
                    drawleg(x, y, leg, :fill)
                end
            end
        end
    end
    finish()
    println("SVG image saved at $(output_path)")
    preview()
end

function create_W_something()
    i = Index(3, "i")
    j = Index(3, "j")
    k = Index(3, "k")
    l = Index(3, "l")
    W = ITensor(Float64, i, j, k, l)

    # W must be symmetric

    # white
    W[i=>1, j=>1, k=>1, l=>1] = 1.5
    # W[i=>1, j=>1, k=>1, l=>1] = 3
    # gray bar
    add_symmetrically(W, (i=>2, j=>1, k=>2, l=>1), 1)
    # gray corner
    add_symmetrically(W, (i=>2, j=>2, k=>1, l=>1), 2)
    # black bar
    add_symmetrically(W, (i=>3, j=>1, k=>3, l=>1), 1)
    # black corner
    add_symmetrically(W, (i=>3, j=>3, k=>1, l=>1), 1.5)
    # branching
    add_symmetrically(W, (i=>2, j=>3, k=>1, l=>3), 2)

    return W
end

function main()

    chi = 30
    maxsize = 30

    width  = 60
    height = 60

    W = create_W_something()
    # run simulation
    println("Starting simulation...")
    start_time = time()
    fixed_ijkl = ctmrg_snapshot(;
        W=W,
        maxsize=maxsize,
        chi=chi,
        width=width,
        height=height
    )
    elapsed = time() - start_time
    println("Execution time: $(round(elapsed, digits=2)) seconds")
    
    # draw snapshot
    draw_something(fixed_ijkl; )
end

main()