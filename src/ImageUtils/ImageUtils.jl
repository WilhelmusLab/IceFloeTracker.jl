"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, imcomplement

include("mask.jl")
include("imcomplement.jl")

end
