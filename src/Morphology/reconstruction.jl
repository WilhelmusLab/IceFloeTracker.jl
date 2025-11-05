import Images: mreconstruct, dilate
import ..skimage: sk_morphology
import ..ImageUtils: to_uint8, imcomplement
import PythonCall: Py, pyconvert

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
    type == "dilation" && return reconstruct(img, se, Val(:dilation), invert)
    type == "erosion" && return reconstruct(img, se, Val(:erosion), invert)
    throw(ArgumentError("Invalid type: $type. Must be 'dilation' or 'erosion'."))
end

function reconstruct(img, se, type::Val{:erosion}, invert::Bool=true)
    morphed =
        sk_morphology.erosion(Py(img).to_numpy(); footprint=Py(collect(se)).to_numpy()) |>
        ((pyimg) -> pyconvert(Array, pyimg)) |>
        to_uint8
    invert && (morphed=imcomplement(to_uint8(morphed)); img=imcomplement(img))
    return pyconvert(
        Array, sk_morphology.reconstruction(Py(morphed).to_numpy(), Py(img).to_numpy())
    )
end

function reconstruct(img, se, type::Val{:dilation}, invert::Bool=true)
    morphed =
        sk_morphology.dilation(Py(img).to_numpy(); footprint=Py(collect(se)).to_numpy()) |>
        ((pyimg) -> pyconvert(Array, pyimg)) |>
        to_uint8
    invert && (morphed=imcomplement(to_uint8(morphed)); img=imcomplement(img))
    return mreconstruct(dilate, morphed, img)
end

"""
    impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}

Use morphological reconstruction to enforce minima on the input image `I` at the positions where the binary mask `BW` is non-zero.

It supports both integer and grayscale images using different implementations for each.
"""
function impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}
    marker = 255 .* BW
    mask = imcomplement(min.(I .+ 1, 255 .- marker))
    reconstructed = pyconvert(
        Array, sk_morphology.reconstruction(Py(marker).to_numpy(), Py(mask).to_numpy())
    )
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

    return 1 .- pyconvert(
        Array,
        sk_morphology.reconstruction(Py(1 .- marker).to_numpy(), Py(1 .- mask).to_numpy()),
    )
end
