"""
    create_landmask(landmask_image, struct_elem, fill_value_lower, fill_value_upper)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix, and use a structuring element to extend a buffer to mask complex coastal features. In the resulting mask, land = 0 and ocean = 1. Returns a named tuple with the dilated and non-dilated landmask.

# Arguments
- `landmask_image`: RGB land mask image from `fetchdata`
- `struct_elem`: structuring element for dilation (optional)
- `fill_value_lower`: fill holes having at least these many pixels (optional)
- `fill_value_upper`: fill holes having at most these many pixels (optional)

"""
function create_landmask(
    landmask_image::T,
    struct_elem::AbstractMatrix{Bool};
    fill_value_lower::Int=0,
    fill_value_upper::Int=2000,
) where {T<:AbstractMatrix}
    landmask_binary = binarize_landmask(landmask_image)
    dilated = dilate(landmask_binary, centered(struct_elem))
    return (
        dilated=ImageMorphology.imfill(.!dilated, (fill_value_lower, fill_value_upper)),
        non_dilated=.!landmask_binary,
    )
end

function create_landmask(landmask_image)
    return create_landmask(landmask_image, make_landmask_se())
end

"""
    binarize_landmask(landmask_image)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix; land = 0, ocean = 1.

# Arguments
- `landmask_image`: RGB land mask image from `fetchdata`
"""
function binarize_landmask(landmask_image::T)::BitMatrix where {T<:AbstractMatrix}
    return ifelse(
        !(typeof(landmask_image) <: AbstractMatrix{Bool}),
        Gray.(landmask_image) .> 0,
        landmask_image,
    )
end

"""
    apply_landmask(input_image, landmask_binary)

Zero out pixels in all channels of the input image using the binary landmask.


# Arguments
- `input_image`: truecolor RGB image
- `landmask_binary`: binary landmask with 1=land, 0=water/ice

"""
function apply_landmask(input_image::AbstractMatrix, landmask_binary::BitMatrix)
    image_masked = landmask_binary .* input_image
    return image_masked
end

# in-place version
function apply_landmask!(input_image::AbstractMatrix, landmask_binary::BitMatrix)
    input_image .= landmask_binary .* input_image
    return nothing
end

"""
    apply_landmask(img, landmask; as_indices::Bool=false)

Apply the landmask to the input image, optionally returning the indices of non-masked (ocean/ice) pixels.

# Arguments
- `img`: input image (e.g., ice mask or RGB image)
- `landmask`: binary landmask (1=ocean/ice, 0=land)
- `as_indices`: if true, return indices of non-masked pixels; otherwise, return masked image
"""
function apply_landmask(img, landmask; as_indices::Bool)
    landmasked = apply_landmask(img, landmask)
    return as_indices ? findall(vec(landmasked)) : landmasked
end
