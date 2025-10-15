import StatsBase: StatsBase
import Images: LinearStretching, adjust_histogram

# TODO: make imadjust more general to work with Gray{Float64} and RGB{Float64} types
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
