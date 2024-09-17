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
