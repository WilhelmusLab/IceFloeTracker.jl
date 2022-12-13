"""
    normalize_image(image_sharpened, image_sharpened_gray, landmask, struct_elem)

Adjusts sharpened land-masked image to highlight ice floe feature.

Does dilation, opening, and landmasking to `image_sharpened`.

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
    struct_elem::Matrix{Bool};
)::Matrix{Gray{Float64}} where {T<:AbstractMatrix{Gray{Float64}}}
    image_dilated = ImageMorphology.dilate(image_sharpened_gray; dims=struct_elem)

    image_opened = ImageMorphology.opening(
        complement.(image_dilated); dims=complement.(image_sharpened)
    )
    return IceFloeTracker.apply_landmask(image_opened, landmask)
end

function normalize_image(
    image_sharpened::Matrix{Float64},
    image_sharpened_gray::AbstractMatrix{Gray{Float64}},
    landmask::BitMatrix,
)::Matrix{Gray{Float64}}
    return normalize_image(
        image_sharpened, image_sharpened_gray, landmask, collect(strel_diamond((5, 5)))
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
        AdaptiveEqualization(;
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
    imsharpen(truecolor_image, lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Sharpen `truecolor_image`.

# Arguments
- `truecolor_image`: input image in truecolor
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
    truecolor_image,
    lambda::Real=0.25,
    kappa::Real=75,
    niters::Int64=3,
    nbins::Int64=255,
    rblocks::Int64=8,
    cblocks::Int64=8,
    clip::Float64=0.8,
    smoothing_param::Int64=10,
    intensity::Float64=2.0,
)::Matrix{Float64}

    image_diffused = diffusion(truecolor_image, lambda, kappa, niters)
    image_diffused_RGB = RGB.(image_diffused)
    masked_view = Float64.(channelview(image_diffused_RGB))

    eq = [
        _adjust_histogram(masked_view[i, :, :], nbins, rblocks, cblocks, clip) for i in 1:3
    ]
    image_equalized = colorview(RGB, eq...)
    image_equalized_gray = Gray.(image_equalized)

    image_smoothed = imfilter(image_equalized_gray, Kernel.gaussian(smoothing_param))

    image_sharpened =
        image_equalized_gray .* (1 + intensity) .+ image_smoothed.* (-intensity)
    image_sharpened = max.(image_sharpened, 0.0)
    return min.(image_sharpened, 1.0)
end

"""
    imsharpen_gray(imgsharpened, landmask)

Apply landmask and return Gray type image in colorview for normalization.
    
"""
function imsharpen_gray(
    imgsharpened::Matrix{Float64}, landmask::AbstractArray{Bool}
)::AbstractMatrix{Gray{Float64}}
    image_sharpened_landmasked = apply_landmask(imgsharpened, landmask)
    return colorview(Gray, image_sharpened_landmasked)
end
