"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, imcomplement, get_brighten_mask, imbrighten, to_uint8

include("mask.jl")
include("imcomplement.jl")
include("brighten.jl")
include("uint8.jl")

end
