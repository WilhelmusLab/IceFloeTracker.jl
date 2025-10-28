"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, get_area_missed, get_optimal_tile_size, get_tile_dims, get_tiles

include("mask.jl")
include("tiling.jl")

end
