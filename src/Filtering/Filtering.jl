module Filtering

export exponential,
    inverse_quadratic,
    SupportedFunctions,
    is_supported,
    SUPPORTED_GRADIENT_FUNCTIONS,
    to_uint8,
    conditional_histeq,
    histeq,
    rgb2gray,
    PeronaMalikDiffusion,
    nonlinear_diffusion

include("gradient_functions.jl")
include("histogram_equalization.jl")
include("imadjust.jl")
include("nonlinear_diffusion.jl")

end
