import Images: Kernel, imfilter, channelview, colorview, base_colorant_type, float64

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
    threshold::Real=0.01,
)
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
    diff_magnitude = if length(size(diff)) > 2
        dropdims(sqrt.(sum(diff .^ 2; dims=1)); dims=1)
    else
        abs.(diff)
    end

    threshold > 0 &&
        (sharpened_image[diff_magnitude .< threshold] .= img[diff_magnitude .< threshold])
    return sharpened_image
end

# method for float matrix with maximum value less than or equal to 1
# TODO: remove after refactor to normed images
function unsharp_mask(img::Matrix{Float64}, smoothing_param, intensity)
    image_sharpened = unsharp_mask(Gray.(img), smoothing_param, intensity, 0)
    return Float64.(image_sharpened)
end

# method for integer matrices
# TODO: this function's `clampmax` argument is unused
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
function unsharp_mask(img::Matrix{Int64}, smoothing_param, intensity, clampmax)
    image_gray = Gray.(img ./ 255)
    image_sharpened = unsharp_mask(image_gray, smoothing_param, intensity, 0)

    return round.(Int, Float64.(image_sharpened) .* 255)
end
