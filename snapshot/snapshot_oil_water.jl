using ITensors, Plots, StatsBase, Luxor
include("snapshot.jl")

function draw_oil_water(fixed_ijkl; cellsize=20, filename="oil_and_water")
    # This function will draw (width-2) by (height-2) diagram

    height, width = size(fixed_ijkl)

    output_path = get_snapshot_output_path(@__DIR__, filename*".svg")

    Drawing((width-2) * cellsize, (height-2) * cellsize, output_path)

    origin(Point(-cellsize, -cellsize)) # origin at top left
    background("white")

    sethue("black")

    legwidth = cellsize / 4

    roundfactor = cellsize / 10

    function drawleg(x, y, n)
        if n == 1
            box(Point(x+cellsize/4+legwidth/4-roundfactor/2, y+cellsize/2), cellsize/2+legwidth/2+roundfactor, legwidth, roundfactor, action = :fill)
        elseif n == 2
            box(Point(x+cellsize/2, y+cellsize/4+legwidth/4-roundfactor/2), legwidth, cellsize/2+legwidth/2+roundfactor, roundfactor, action = :fill)
        elseif n == 3
            box(Point(x+cellsize/4*3-legwidth/4+roundfactor/2, y+cellsize/2), cellsize/2+legwidth/2+roundfactor, legwidth, roundfactor, action = :fill)
        elseif n == 4
            box(Point(x+cellsize/2, y+cellsize/4*3-legwidth/4+roundfactor/2), legwidth, cellsize/2+legwidth/2+roundfactor, roundfactor, action = :fill)
        end
    end

    for row in 2:height-1
        for col in 2:width-1
            x = (col - 1) * cellsize
            y = (row - 1) * cellsize
            ijkl = fixed_ijkl[row, col]

            if ijkl == (2,2,2,2)
                box(Point(x+cellsize/2, y+cellsize/2), cellsize, cellsize, roundfactor, action = :fill)

                around = [x==(2,2,2,2) for x in [fixed_ijkl[row, col-1], fixed_ijkl[row-1, col], fixed_ijkl[row, col+1], fixed_ijkl[row+1, col]]]

                if around[1]
                    box(Point(x, y), Point(x+roundfactor, y+cellsize), action = :fill)
                end
                if around[2]
                    box(Point(x, y), Point(x+cellsize, y+roundfactor), action = :fill)
                end
                if around[3]
                    box(Point(x+cellsize-roundfactor, y), Point(x+cellsize, y+cellsize), action = :fill)
                end
                if around[4]
                    box(Point(x, y+cellsize-roundfactor), Point(x+cellsize, y+cellsize), action = :fill)
                end
            else
                for leg in 1:4
                    if ijkl[leg] == 2
                        drawleg(x, y, leg)
                    end
                end
            end
        end
    end
    finish()
    println("SVG image saved at $(output_path)")
    preview()
end

function CUI_draw_oil_water(fixed_ijkl)
    # CUI output
    
    ijkl_unicode = Dict(
        (1,1,1,1) => " ",
        (1,1,1,2) => "╷",
        (1,1,2,1) => "╶",
        (1,2,1,1) => "╵",
        (2,1,1,1) => "╴",

        (1,1,2,2) => "┌",
        (1,2,2,1) => "└",
        (2,2,1,1) => "┘",
        (2,1,1,2) => "┐",

        (2,2,2,2) => "█",
    )
    height, width = size(fixed_ijkl)
    for row in 1:height
        println()
        for col in 1:width
            ijkl = fixed_ijkl[row, col]
            print(ijkl_unicode[ijkl])
        end
    end
    println()
end

function create_W_oil_water(x, a)
    i = Index(2, "i")
    j = Index(2, "j")
    k = Index(2, "k")
    l = Index(2, "l")
    W = ITensor(Float64, i, j, k, l)

    # W must be symmetric
    W[i=>1, j=>1, k=>1, l=>1] = 1.0
    W[i=>2, j=>1, k=>1, l=>1] = x
    W[i=>1, j=>2, k=>1, l=>1] = x
    W[i=>1, j=>1, k=>2, l=>1] = x
    W[i=>1, j=>1, k=>1, l=>2] = x
    W[i=>2, j=>2, k=>1, l=>1] = a * x^2
    W[i=>1, j=>2, k=>2, l=>1] = a * x^2
    W[i=>1, j=>1, k=>2, l=>2] = a * x^2
    W[i=>2, j=>1, k=>1, l=>2] = a * x^2
    W[i=>2, j=>2, k=>2, l=>2] = x^4

    return W
end

function main()
    # Define constants
    a = 0.3599
    x = 1.31438

    chi = 30
    maxsize = 30

    width  = 100
    height = 100

    # run simulation
    println("Starting simulation...")
    start_time = time()
    fixed_ijkl = ctmrg_snapshot(;
        W=create_W_oil_water(x, a),
        maxsize=maxsize,
        chi=chi,
        width=width,
        height=height
    )
    elapsed = time() - start_time
    println("Execution time: $(round(elapsed, digits=2)) seconds")

    # draw snapshot
    # CUI output
    CUI_draw_oil_water(fixed_ijkl)
    # SVG output
    draw_oil_water(fixed_ijkl; filename="oilwater_a$(a)_x$(x)_chi$(chi)_N$(maxsize)_W$(width)_H$(height)")
end

main()