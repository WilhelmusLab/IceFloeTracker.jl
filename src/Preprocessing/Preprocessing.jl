module Preprocessing

abstract type IceFloePreprocessingAlgorithm end

struct LopezAcostaPreprocessing2019 <: IceFloePreprocessingAlgorithm
    param1
    param2
end

struct WatkinsPreprocessing2025 <: IceFloePreprocessingAlgorithm
    param1
    param2
end

function WatkinsPreprocessing2025(; landmask_structuring_element=make_landmask_se())
    return WatkinsPreprocessing2025(landmask_structuring_element)
end

"""
Produce a grayscale image with sharpened floe edges and smoothed floe interiors, following Watkins et al. 2025.
"""
function (p::WatkinsPreprocessing2025)
    (truecolor_image,
     falsecolor_image,
     cloud_mask,
     land_mask,
     tiles)
    # image diffusion
    # equalization
    # sharpening (local)
    # morphological residue
    # gamma adjustment
    # sharpening (regional)
end


function LopezAcostaPreprocessing2019(; landmask_structuring_element=make_landmask_se())
    return LopezAcosta2019(landmask_structuring_element)
end

function (p::LopezAcostaPreprocessing2019)()
    # initial algorithm steps
end






export LopezAcostaPreprocessing2019
export WatkinsPreprocessing2025

end
