"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, imcomplement, to_uint8

include("mask.jl")
include("imcomplement.jl")
include("uint8.jl")

end
