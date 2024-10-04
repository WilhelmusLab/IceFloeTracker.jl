"""
    normalize_image(image_sharpened, image_sharpened_gray, landmask, struct_elem;)

Adjusts sharpened land-masked image to highlight ice floe features.

Does reconstruction and landmasking to `image_sharpened`.

# Arguments
- `image_sharpened`: sharpened image (output of `imsharpen`)
- `image_sharpened_gray`: grayscale, landmasked sharpened image (output of `imsharpen_gray(image_sharpened)`)
- `landmask`: landmask for region of interest
- `struct_elem`: structuring element for dilation

"""
function normalize_image(
    image_sharpened::Matrix{Float64},
    image_sharpened_gray::T,
    landmask::BitMatrix,
    struct_elem::ImageMorphology.MorphologySEArray{2};
)::Matrix{Gray{Float64}} where {T<:AbstractMatrix{Gray{Float64}}}
    image_dilated = MorphSE.dilate(image_sharpened_gray, struct_elem)

    image_reconstructed = MorphSE.mreconstruct(
        MorphSE.dilate, complement.(image_dilated), complement.(image_sharpened)
    )
    return IceFloeTracker.apply_landmask(image_reconstructed, landmask)
end

function normalize_image(
    image_sharpened::Matrix{Float64},
    image_sharpened_gray::Matrix{Gray{Float64}},
    landmask::BitMatrix,
)::Matrix{Gray{Float64}}
    return normalize_image(
        image_sharpened, image_sharpened_gray, landmask, strel_diamond((5, 5))
    )
end

"""
    _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)

Perform adaptive histogram equalization to a masked image. To be invoked within `imsharpen`.

# Arguments
- `masked_view`: input image in truecolor
See `imsharpen` for a description of the remaining arguments

"""
function _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)
    return adjust_histogram(
        masked_view,
        ImageContrastAdjustment.AdaptiveEqualization(;
            nbins=nbins,
            rblocks=rblocks,
            cblocks=cblocks,
            minval=minimum(masked_view),
            maxval=maximum(masked_view),
            clip=clip,
        ),
    )
end

"""
    imsharpen(truecolor_image, landmask_no_dilate, lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Sharpen `truecolor_image`.

# Arguments
- `truecolor_image`: input image in truecolor
- `landmask_no_dilate`: landmask for region of interest
- `lambda`: speed of diffusion (0–0.25)
- `kappa`: conduction coefficient for diffusion (25–100)
- `niters`: number of iterations of diffusion
- `nbins`: number of bins during histogram equalization
- `rblocks`: number of row blocks to divide input image during equalization
- `cblocks`: number of column blocks to divide input image during equalization
- `clip`: Thresholds for clipping histogram bins (0–1); values closer to one minimize contrast enhancement, values closer to zero maximize contrast enhancement
- `smoothing_param`: pixel radius for gaussian blurring (1–10)
- `intensity`: amount of sharpening to perform
"""
function imsharpen(
    truecolor_image::Matrix{RGB{Float64}},
    landmask_no_dilate::BitMatrix,
    lambda::Real=0.1,
    kappa::Real=75,
    niters::Int64=3,
    nbins::Int64=255,
    rblocks::Int64=10, # matlab default is 8 CP
    cblocks::Int64=10, # matlab default is 8 CP
    clip::Float64=0.86, # matlab default is 0.01 CP
    smoothing_param::Int64=10,
    intensity::Float64=2.0,
)::Matrix{Float64}
    input_image = IceFloeTracker.apply_landmask(truecolor_image, landmask_no_dilate)

    input_image .= IceFloeTracker.diffusion(input_image, lambda, kappa, niters)

    masked_view = Float64.(channelview(input_image))

    eq = [
        _adjust_histogram(@view(masked_view[i, :, :]), nbins, rblocks, cblocks, clip) for
        i in 1:3
    ]

    image_equalized = colorview(RGB, eq...)

    image_equalized_gray = Gray.(image_equalized)

    return unsharp_mask(image_equalized_gray, smoothing_param, intensity)
end

"""
    unsharp_mask(image_gray, smoothing_param, intensity, clampmax)

Apply unsharp masking on (equalized) grayscale ([0, `clampmax`]) image to enhance its sharpness.

# Arguments
- `image_gray`: The input grayscale image, typically already equalized.
- `smoothing_param::Int`: The pixel radius for Gaussian blurring (typically between 1 and 10).
- `intensity`: The amount of sharpening to apply. Higher values result in more pronounced sharpening.
- `clampmax`: upper limit of intensity values in the returned image.`
# Returns
The sharpened grayscale image with values clipped between 0 and `clapmax`.
"""
function unsharp_mask(image_gray, smoothing_param, intensity, clampmax)
    image_smoothed = imfilter(image_gray, Kernel.gaussian(smoothing_param))
    clamp!(image_smoothed, 0.0, clampmax)
    image_sharpened = image_gray * (1 + intensity) .- image_smoothed * intensity
    clamp!(image_sharpened, 0.0, clampmax)
    return round.(Int, image_sharpened)
end

# For old workflow in final2020.m
"""
    (Deprecated)
    unsharp_mask(image_gray, smoothing_param, intensity)

Apply unsharp masking on (equalized) grayscale image to enhance its sharpness.

Does not perform clamping after the smoothing step. Kept for legacy tests of IceFloeTracker.jl.

# Arguments
- `image_gray`: The input grayscale image, typically already equalized.
- `smoothing_param::Int`: The pixel radius for Gaussian blurring (typically between 1 and 10).
- `intensity`: The amount of sharpening to apply. Higher values result in more pronounced sharpening.
# Returns
The sharpened grayscale image with values clipped between 0 and `clapmax`.
"""
function unsharp_mask(image_gray, smoothing_param, intensity)
    image_smoothed = imfilter(image_gray, Kernel.gaussian(smoothing_param))
    image_sharpened = image_gray * (1 + intensity) .- image_smoothed * intensity
    return clamp.(image_sharpened, 0.0, 1.0)
end

"""
    imsharpen_gray(imgsharpened, landmask)

Apply landmask and return Gray type image in colorview for normalization.

"""
function imsharpen_gray(
    imgsharpened::Matrix{Float64}, landmask::AbstractArray{Bool}
)::Matrix{Gray{Float64}}
    image_sharpened_landmasked = apply_landmask(imgsharpened, landmask)
    return colorview(Gray, image_sharpened_landmasked)
end


