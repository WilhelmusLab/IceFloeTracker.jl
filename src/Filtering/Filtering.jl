module Filtering

export exponential,
    inverse_quadratic,
    SupportedFunctions,
    is_supported,
    SUPPORTED_GRADIENT_FUNCTIONS,
    conditional_histeq,
    histeq,
    rgb2gray,
    nonlinear_diffusion,
    AbstractDiffusionAlgorithm,
    PeronaMalikDiffusion,
    anisotropic_diffusion_3D,
    anisotropic_diffusion_2D,
    adapthisteq,
    imadjust,
    unsharp_mask

include("gradient_functions.jl")
include("histogram_equalization.jl")
include("imadjust.jl")
include("nonlinear_diffusion.jl")
include("unsharp_mask.jl")

end
