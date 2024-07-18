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
    rblocks::Int64=10,
    cblocks::Int64=10,
    clip::Float64=0.86,
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

    image_smoothed = imfilter(image_equalized_gray, Kernel.gaussian(smoothing_param))

    image_sharpened =
        image_equalized_gray .* (1 + intensity) .+ image_smoothed .* (-intensity)
    return min.(max.(image_sharpened, 0.0), 1.0)
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
