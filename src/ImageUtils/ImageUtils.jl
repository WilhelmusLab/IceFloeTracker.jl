# IceFloeTracker.jl exports:
export masker,
    apply_mask,
    get_area_missed,
    get_brighten_mask,
    get_optimal_tile_size,
    get_tile_dims,
    get_tiles,
    imbrighten,
    imcomplement,
    to_uint8

"""
Module for general image utilities.
"""
module ImageUtils

export masker,
    apply_mask,
    get_area_missed,
    get_brighten_mask,
    get_optimal_tile_size,
    get_tile_dims,
    get_tiles,
    imbrighten,
    imcomplement,
    to_uint8

include("brighten.jl")
include("imcomplement.jl")
include("mask.jl")
include("tiling.jl")
include("uint8.jl")

end
