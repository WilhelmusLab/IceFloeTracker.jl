import Images: mreconstruct, dilate, erode
import ..ImageUtils: to_uint8, imcomplement

"""
    reconstruct(img, se, type, invert)

Perform closing/opening by reconstruction on `img`.

# Arguments
- `img::AbstractArray`: The input image.
- `se::AbstractArray`: The structuring element.
- `type::String`: The type of morphological operation to perform. Must be either "dilation" (close by reconstruction) or `"erosion"` (open by reconstruction).
- `invert::Bool=true`: Invert marker and mask before reconstruction.
"""
function reconstruct(img, se, type, invert::Bool=true)
    !(type == "dilation" || type == "erosion") &&
        throw(ArgumentError("Invalid type: $type. Must be 'dilation' or 'erosion'."))

    type == "dilation" && (morphed = to_uint8(dilate(img, collect(se))))
    type == "erosion" && (morphed = to_uint8(erode(img, collect(se))))

    invert && (morphed = imcomplement(to_uint8(morphed)); img = imcomplement(img))

    type == "dilation" && return mreconstruct(dilate, morphed, img)

    return mreconstruct(dilate, morphed, img)
end

"""
    impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}

Use morphological reconstruction to enforce minima on the input image `I` at the positions where the binary mask `BW` is non-zero.
"""
function impose_minima(
    I::AbstractArray{T}, BW::AbstractMatrix{Bool}
) where {T<:AbstractFloat}
    # compute shift
    a, b = extrema(I)
    rng = b - a
    h = rng == 0 ? 0.1 : rng / 1000

    marker = -Inf * BW .+ (Inf * .!BW)
    mask = min.(I .+ h, marker)
    return 1 .- mreconstruct(dilate, 1 .- marker, 1 .- mask)
end
