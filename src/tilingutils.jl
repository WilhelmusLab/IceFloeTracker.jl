function getfit(dims::Tuple{Int,Int}, l::Int)::Tuple{Int,Int}
    return dims .÷ l
end

function get_area_missed(l::Int, dims::Tuple{Int,Int}, area::Int)::Float64
    return 1 - prod(getfit(dims, l)) * l^2 / area
end

"""
    get_optimal_tile_size(l0::Int, dims::Tuple{Int,Int}) -> Int

Calculate the optimal tile size in the range [l0-1, l0+1] for the given size `l0` and image dimensions `dims`.

# Description
This function computes the optimal tile size for tiling an area with given dimensions. It ensures that the initial tile size `l0` is at least 2 and not larger than any of the given dimensions. The function evaluates candidate tile sizes and selects the one that minimizes the area missed during tiling. In case of a tie, it prefers the larger tile size.

# Example
```
julia> get_optimal_tile_size(3, (10, 7))
2
```
"""
function get_optimal_tile_size(l0::Int, dims::Tuple{Int,Int})::Int
    l0 < 2 && error("l0 must be at least 2")
    any(l0 .> dims) && error("l0 = $l0 is too large for the given dimensions $dims")

    area = prod(dims)
    minimal_shift = l0 == 2 ? 0 : 1
    candidates = [l0 + i for i in -minimal_shift:1]

    minl, M = 0, Inf
    for l in candidates
        missedarea = get_area_missed(l, dims, area)
        if missedarea <= M # prefer larger l in case of tie
            M, minl = missedarea, l
        end
    end
    return minl
end

"""
    get_tile_meta(tile)

Extracts metadata from a given tile.

# Arguments
- `tile`: A collection of tuples, where each tuple represents a coordinate pair.

# Returns
- A tuple `(a, b, c, d)` where:
  - `a`: The first element of the first tuple in `tile`.
  - `b`: The last element of the first tuple in `tile`.
  - `c`: The first element of the last tuple in `tile`.
  - `d`: The last element of the last tuple in `tile`.
"""
function get_tile_meta(tile)
    a, c = first.(tile)
    b, d = last.(tile)
    return [a, b, c, d]
end


function bump_tile(tile, dims)
    extrarows, extracols = dims
    a, b, c, d = get_tile_meta(tile)
    b += extrarows
    d += extracols
    return (a:b, c:d)
end

function get_tile_dims(tile)
    a, b, c, d = get_tile_meta(tile)
    width, height = d - c + 1, b - a + 1
    return (width, height)
end

"""
    adjust_edge_tiles(tiles, bumpby)

Adjusts the edge tiles of a tiling by bumping them with the given dimensions.

The algorithm works in two steps: first fold the right edge tiles, then fold the bottom edge tiles.

# Arguments
- `tiles`: A collection of tiles.
- `bumpby`: A tuple `(extrarows, extracols)` representing the dimensions to bump the tiles.
"""
function adjust_edge_tiles(tiles, bumpby=nothing)
    if isnothing(bumpby)
        bumpby = get_tile_dims(tiles[end])
    end

    _, l = get_tile_dims(tiles[1])
    shift_height, shift_width = 0, 0

    # Fold right edge tiles if leftover width is less than half of the standard width
    if last(bumpby) < l ÷ 2
        tiles_right_edge = @view tiles[:, end-1]
        tiles_right_edge .= bump_tile.(tiles_right_edge, Ref((0, last(bumpby))))
        shift_height += 1
    end

    if first(bumpby) < l ÷ 2
        tiles_bottom_edge = @view tiles[end-1, :]
        tiles_bottom_edge .= bump_tile.(tiles_bottom_edge, Ref((first(bumpby), 0)))
        shift_width += 1
    end
    return tiles[1:end-shift_width, 1:end-shift_height]
end