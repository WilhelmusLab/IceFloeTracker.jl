import Images: mreconstruct, dilate
import ..skimage: sk_morphology
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

    type == "dilation" &&
        (morphed = to_uint8(sk_morphology.dilation(img; footprint=collect(se))))
    type == "erosion" &&
        (morphed = to_uint8(sk_morphology.erosion(img; footprint=collect(se))))

    invert && (morphed = imcomplement(to_uint8(morphed)); img = imcomplement(img))

    type == "dilation" && return mreconstruct(dilate, morphed, img)

    return sk_morphology.reconstruction(morphed, img)
end

"""
    impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}

Use morphological reconstruction to enforce minima on the input image `I` at the positions where the binary mask `BW` is non-zero.

It supports both integer and grayscale images using different implementations for each.
"""
function impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}
    marker = 255 .* BW
    mask = imcomplement(min.(I .+ 1, 255 .- marker))
    reconstructed = sk_morphology.reconstruction(marker, mask)
    return imcomplement(Int.(reconstructed))
end

function impose_minima(
    I::AbstractArray{T}, BW::AbstractMatrix{Bool}
) where {T<:AbstractFloat}
    # compute shift
    a, b = extrema(I)
    rng = b - a
    h = rng == 0 ? 0.1 : rng / 1000

    marker = -Inf * BW .+ (Inf * .!BW)
    mask = min.(I .+ h, marker)

    return 1 .- sk_morphology.reconstruction(1 .- marker, 1 .- mask)
end
