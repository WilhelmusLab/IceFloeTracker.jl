"""
Module for general image utilities.
"""
module ImageUtils

export masker,
    apply_mask,
    apply_to_channels,
    get_area_missed,
    get_brighten_mask,
    get_optimal_tile_size,
    get_tile_dims,
    get_tiles,
    imbrighten,
    imcomplement,
    to_uint8,
    get_tile_meta,
    bump_tile

include("brighten.jl")
include("imcomplement.jl")
include("mask.jl")
include("tiling.jl")
include("uint8.jl")
include("misc_image_utils.jl")

end
