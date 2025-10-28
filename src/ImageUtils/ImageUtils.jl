"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, imcomplement, to_uint8

include("imcomplement.jl")
include("mask.jl")
include("uint8.jl")

end
