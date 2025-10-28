"""
Module for general image utilities.
"""
module ImageUtils

export masker, apply_mask, get_brighten_mask, imbrighten

include("mask.jl")
include("brighten.jl")

end
