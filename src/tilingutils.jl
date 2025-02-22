
"""
    getfit(dims::Tuple{Int,Int}, side_length::Int)::Tuple{Int,Int}

Calculate how many tiles of a given side length fit into the given dimensions.

# Arguments
- `dims::Tuple{Int,Int}`: A tuple representing the dimensions (width, height).
- `side_length::Int`: The side length of the tile.

# Returns
- `Tuple{Int,Int}`: A tuple representing the number of tiles that fit along each dimension.

# Examples
```
julia> getfit((10, 20), 5)
(2, 4)

julia> getfit((15, 25), 5)
(3, 5)
"""
function getfit(dims::Tuple{Int,Int}, side_length::Int)::Tuple{Int,Int}
    return dims .÷ side_length
end


"""
    get_area_missed(side_length::Int, dims::Tuple{Int,Int})::Float64

Calculate the proportion of the area that is not covered by tiles of a given side length.

# Arguments
- `side_length::Int`: The side length of the tile.
- `dims::Tuple{Int,Int}`: A tuple representing the dimensions (width, height).

# Returns
- `Float64`: The proportion of the area that is not covered by the tiles.

# Examples
```
julia> get_area_missed(5, (10, 20))
0.0

julia> get_area_missed(7, (10, 20))
0.51
"""
function get_area_missed(side_length::Int, dims::Tuple{Int,Int})::Float64
    area = prod(dims)
    return 1 - prod(getfit(dims, side_length)) * side_length^2 / area
end


"""
    get_optimal_tile_size(l0::Int, dims::Tuple{Int,Int}) -> Int

Calculate the optimal tile size in the range [l0-1, l0+1] for the given size `l0` and image dimensions `dims`.

# Description
This function computes the optimal tile size for tiling an area with given dimensions. It ensures that the initial tile size `l0` is at least 2 and not larger than any of the given dimensions. The function evaluates candidate tile sizes and selects the one that minimizes the area missed by its corresponding tiling. In case of a tie, it prefers the larger tile size.

# Example
```
julia> get_optimal_tile_size(3, (10, 7))
2
```
"""
function get_optimal_tile_size(l0::Int, dims::Tuple{Int,Int})::Int
    l0 < 2 && error("l0 must be at least 2")
    any(l0 .> dims) && error("l0 = $l0 is too large for the given dimensions $dims")

    minimal_shift = l0 == 2 ? 0 : 1
    candidates = [l0 + i for i in -minimal_shift:1]

    minl, M = 0, Inf
    for side_length in candidates
        missedarea = get_area_missed(side_length, dims)
        if missedarea <= M # prefer larger side_length in case of tie
            M, minl = missedarea, side_length
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

"""
    bump_tile(tile::Tuple{UnitRange{Int64}, UnitRange{Int64}}, dims::Tuple{Int,Int})::Tuple{UnitRange{Int}, UnitRange{Int}}

Adjust the tile dimensions by adding extra rows and columns.

# Arguments
- `tile::Tuple{Int,Int,Int,Int}`: A tuple representing the tile dimensions (a, b, c, d).
- `dims::Tuple{Int,Int}`: A tuple representing the extra rows and columns to add (extrarows, extracols).

# Returns
- `Tuple{UnitRange{Int}, UnitRange{Int}}`: A tuple of ranges representing the new tile dimensions.

# Examples
```julia
julia> bump_tile((1:3, 1:4), (1, 1))
(1:4, 1:5)
"""
function bump_tile(tile::Tuple{UnitRange{S},UnitRange{S}}, dims::Tuple{S,S}) where {S<:Int}
    extrarows, extracols = dims
    a, b, c, d = get_tile_meta(tile)
    b += extrarows
    d += extracols
    return (a:b, c:d)
end

"""
    get_tile_dims(tile)

Calculate the dimensions of a tile.

# Arguments
- `tile::Tuple{UnitRange{Int},UnitRange{Int}}`: A tuple representing the tile dimensions.

# Returns
- `Tuple{Int,Int}`: A tuple representing the width and height of the tile.

# Examples
```julia
julia> get_tile_dims((1:3, 1:4))
(4, 3)
"""
function get_tile_dims(tile)
    a, b, c, d = get_tile_meta(tile)
    width, height = d - c + 1, b - a + 1
    return (width, height)
end


"""
    get_tiles(array, t::Tuple{Int,Int})

Generate a collection of tiles from an array.

The function adjusts the bottom and right edges of the tile matrix if they are smaller than half the tile sizes in `t`.
"""
function get_tiles(array, t::Tuple{T,T}) where T<:Union{Int,Int64}
    a, b = t
    tiles = TileIterator(axes(array), (a, b)) |> collect
    _a, _b = size(array)

    bottombump = mod(_a, a)
    rightbump = mod(_b, b)

    if bottombump == 0 && rightbump == 0
        return tiles
    end

    crop_height, crop_width = 0, 0

    # Adjust bottom edge if necessary
    if bottombump <= a ÷ 2
        bottom_edge = tiles[end-1, :]
        tiles[end-1, :] .= bump_tile.(bottom_edge, Ref((bottombump, 0)))
        crop_height += 1
    end

    # Adjust right edge if necessary
    if rightbump <= b ÷ 2
        right_edge = tiles[:, end-1]
        tiles[:, end-1] .= bump_tile.(right_edge, Ref((0, rightbump)))
        crop_width += 1
    end

    return tiles[1:end-crop_height, 1:end-crop_width]
end

"""
    get_tiles(array, side_length)

Generate a collection of tiles from an array.

Unlike `TileIterator`, the function adjusts the bottom and right edges of the tile matrix if they are smaller than half the tile size `side_length`.
"""
function get_tiles(array, side_length::Int)
    return get_tiles(array, (side_length, side_length))
end

function get_tiles(array; rblocks, cblocks)
    rtile, ctile = size(array)
    tile_size = (rtile ÷ rblocks, ctile ÷ cblocks)
    return TileIterator(axes(array), tile_size)
end
