# Helper functions

"""
    loadimg(; dir::String, fname::String)

Load an image from `dir` with filename `fname` into a matrix of `Float64` values. Returns the loaded image.
"""
function loadimg(; dir::String, fname::String)
    return (x -> float64.(x))(load(joinpath(dir, fname)))
end

"""
    add_padding(img, style)

Extrapolate the image `img` according to the `style` specifications type. Returns the extrapolated image.

# Arguments
- `img`: Image to be padded.
- `style`: A supported type (such as `Pad` or `Fill`) representing the extrapolation style. See the relevant [documentation](https://juliaimages.org/latest/function_reference/#ImageFiltering) for details.

See also [`remove_padding`](@ref)
"""
function add_padding(img, style::Union{Pad,Fill})::Matrix
    return collect(Images.padarray(img, style))
end

"""
    remove_padding(paddedimg, border_spec)

Removes padding from the boundary of padded image `paddedimg` according to the border specification `border_spec` type. Returns the cropped image.

# Arguments
- `paddedimg`: Pre-padded image.
- `border_spec`: Type representing the style of padding (such as `Pad` or `Fill`) with which `paddedimg` is assumend to be pre-padded. Example: `Pad((1,2), (3,4))` specifies 1 row on the top, 2 columns on the left, 3 rows on the bottom, and 4 columns on the right boundary.

See also [`add_padding`](@ref)
"""
function remove_padding(paddedimg, border_spec::Union{Pad,Fill})::Matrix
    top, left = border_spec.lo
    bottom, right = border_spec.hi
    return paddedimg[(top + 1):(end - bottom), (left + 1):(end - right)]
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
    return IceFloeTracker.imcomplement(Int.(reconstructed))
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

"""
    bwdist(bwimg)

Distance transform for binary image `bwdist`.
"""
function bwdist(bwimg::AbstractArray{Bool})::AbstractArray{Float64}
    return Images.distance_transform(Images.feature_transform(bwimg))
end
