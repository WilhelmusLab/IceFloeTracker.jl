"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, create_landmask, apply_landmask, apply_landmask!

include("landmask.jl")
include("mask.jl")

end
