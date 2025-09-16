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
    struct_elem;
)::Matrix{Gray{Float64}} where {T<:AbstractMatrix{Gray{Float64}}}
    image_dilated = dilate(image_sharpened_gray, struct_elem)

    image_reconstructed = mreconstruct(
        dilate, complement.(image_dilated), complement.(image_sharpened)
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
    kappa::Real=0.1,
    niters::Int64=3,
    nbins::Int64=255,
    rblocks::Int64=10, # matlab default is 8 CP
    cblocks::Int64=10, # matlab default is 8 CP
    clip::Float64=0.86, # matlab default is 0.01 CP
    smoothing_param::Int64=10,
    intensity::Float64=2.0,
)::Matrix{Float64}
    input_image = IceFloeTracker.apply_landmask(truecolor_image, landmask_no_dilate)

    input_image .= IceFloeTracker.nonlinear_diffusion(input_image, lambda, kappa, niters)

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


# TODO: Remove function, replace with direct use of landmask and colorview.
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

# TODO: Remove once the workflow is all normed images
function adjustgamma(img, gamma=1.5, asuint8=true)
    if maximum(img) > 1
        img = img ./ 255
    end

    adjusted = adjust_histogram(img, GammaCorrection(gamma))

    if asuint8
        adjusted = Int.(round.(adjusted * 255, RoundNearestTiesAway))
    end

    return adjusted
end

# TODO: Remove function
function imbinarize(img)
    f = AdaptiveThreshold(img) # infer the best `window_size` using `img`
    return binarize(img, f)
end


# TODO: Move to a new module for image filters
"""
    unsharp_mask(img, radius, amount, threshold)

    Enhance image sharpness by weighted differencing of the image and a Gaussian blurred image.
    If ``B`` is the blurred version of image ``I``, then an unsharp mask sharpened image is obtained by
    ``S = I + (I - B)*A``
    The amount of sharpening is determined by the factor A. An option threshold can be supplied such
    that the sharpening is only applied where ``I - B`` is greater than some factor.

    # Arguments
    img: input image
    radius: standard deviation of the Gaussian blur
    amount: multiplicative factor
    threshold: minimum difference for applying the sharpening

    # Returns
    Sharpened image
"""
function unsharp_mask(
    img::AbstractArray{<:Union{AbstractRGB,TransparentRGB,AbstractGray}},
    radius::Real=3,
    amount::Real=0.5,
    threshold::Real=0.01)

    image_float = float64.(img)
    image_smoothed = imfilter(image_float, Kernel.gaussian(radius))

    cv_image = channelview(image_float)
    cv_smooth = channelview(image_smoothed)    
    diff = cv_image .- cv_smooth
    cv_sharp = cv_image .+ diff .* amount
    clamp!(cv_sharp, 0, 1)

    # Convert back into an image of the original type
    sharpened_image = colorview(base_colorant_type(eltype(img)), cv_sharp)
    sharpened_image = convert.(eltype(img), sharpened_image)

    # Optionally use only where sharpening is larger than threshold
    diff_magnitude = length(size(diff)) > 2 ? dropdims(sqrt.(sum(diff.^2, dims=1)); dims=1) : abs.(diff)

    threshold > 0 && (sharpened_image[diff_magnitude .< threshold] .= img[diff_magnitude .< threshold])
    return sharpened_image
end

# method for float matrix with maximum value less than or equal to 1
# TODO: remove after refactor to normed images
function unsharp_mask(img::Matrix{Float64}, smoothing_param, intensity)
    image_sharpened = unsharp_mask(Gray.(img), smoothing_param, intensity, 0)
    return Float64.(image_sharpened)
end

# method for integer matrices
function unsharp_mask(img::Matrix{Int64}, smoothing_param, intensity, clampmax)
    image_gray = Gray.(img ./ 255)
    image_sharpened = unsharp_mask(image_gray, smoothing_param, intensity, 0)

    return round.(Int, Float64.(image_sharpened) .* 255)
end
