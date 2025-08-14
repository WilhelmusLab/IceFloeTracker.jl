# TODO: make imadjust more general to work with Gray{Float64} and RGB{Float64} types
# dmw: This is a good example of a function where we should only very lightly touch the ImageContrastAdjustment function.
# We don't need to do any image conversion: it should match the input and output types. We can restrict it to grayscale images
# or single channels. The wrapper is just to add percentile-based max and min for the linear stretching.
# Perhaps the new function can be PercentileLinearStretching?
"""
    imadjust(img; low, high)

Adjust the contrast of an image using linear stretching. The image is normalized to [0, 1] and then stretched to the range [low, high].

# Arguments
- `img`: The input image.
- `low`: The lower bound of the stretched image. Default is 0.01.
- `high`: The upper bound of the stretched image. Default is 0.99.

# Returns

The contrast-adjusted image in the range [0, 255].
"""
function imadjust(
    img::AbstractArray{<:Integer}; low::T=0.01, high::T=0.99
)::Matrix{Int} where {T<:AbstractFloat}
    img = img ./ 255
    imgflat = vec(img)
    plow = StatsBase.percentile(imgflat, low * 100)
    phigh = StatsBase.percentile(imgflat, high * 100)

    f = LinearStretching((plow, phigh) => (0.0, 1.0))

    return to_uint8(adjust_histogram(img, f) * 255)
end
