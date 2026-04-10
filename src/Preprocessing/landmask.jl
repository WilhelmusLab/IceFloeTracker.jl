import ..Morphology: _generate_se!
import Images: Gray, dilate, imfill
import OffsetArrays: centered

"""
    make_landmask_se()

Create a non-regular octagonal structuring element for dilating the landmask. 
This structuring element matches a polygonal ``disk'' element used in the MATLAB IFT prototype.
"""
function make_landmask_se()
    se = [sum(c.I) <= 29 for c in CartesianIndices((99, 99))]
    _generate_se!(se)
    return centered(se)
end

"""
    create_landmask(landmask_image, struct_elem, fill_value_lower, fill_value_upper)

Convert a land mask image to a 1-channel binary matrix, and use a structuring element to extend a buffer to mask complex coastal features, and fill holes in the dilated image. In the resulting mask, land = 0 and ocean = 1. Returns a named tuple with the dilated and non-dilated landmask.

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
    land_mask = binarize_landmask(landmask_image)
    coastal_buffer_mask = create_coastal_buffer_mask(
        land_mask,
        struct_elem;
        fill_min_pixels=fill_value_lower,
        fill_max_pixels=fill_value_upper,
    )

    return (dilated=coastal_buffer_mask, non_dilated=land_mask)
end

function create_landmask(landmask_image; strel=make_landmask_se())
    return create_landmask(landmask_image, strel)
end

"""
    create_coastal_buffer_mask(landmask_binary; fill_min_pixels, fill_max_pixels)
    create_coastal_buffer_mask(landmask_binary, structuring_element; fill_min_pixels, fill_max_pixels)


Dilate the binary landmask using the provided structuring element, and fill holes in the dilated image. 
In the input landmask, land = 1 and ocean = 0. 
In the resulting mask, land and coastal buffer = 1, ocean = 0.

The dilation will create a buffer around the land, which can help mask complex coastal features. 
The hole filling step can help fill in small gaps in the dilated mask that may occur due to the dilation process.
"""
function create_coastal_buffer_mask(
    landmask::AbstractMatrix{Bool},
    structuring_element::AbstractMatrix{Bool}=make_landmask_se();
    fill_min_pixels::Int=0,
    fill_max_pixels::Int=2000,
)::Matrix{Bool}
    centered_structuring_element = centered(structuring_element)
    mask_unfilled = dilate(landmask, centered_structuring_element)
    mask_filled = .!imfill(.!mask_unfilled, (fill_min_pixels, fill_max_pixels))
    return mask_filled
end

"""
    binarize_landmask(landmask_image)

Convert a 3-channel RGB or 1-channel Gray land mask image to a 1-channel binary matrix with land = 1, ocean = 0.
Assumes that the input image is 0 over the ocean and some shade over land; the tol argument lets a higher threshold
for land pixels be chosen.

# Arguments
- `landmask_image`: land mask image, e.g. from NASA Worldview
- `tol` (Optional): Values in the image larger than `tol` are considered land.
"""
function binarize_landmask(landmask_image; tol=0.1)::BitMatrix
    return Gray.(landmask_image) .> tol
end

"""
    apply_landmask(input_image, landmask_binary)
    apply_landmask!(input_image, landmask_binary)

Zero out pixels in all channels of the input image using the binary landmask.

# Arguments
- `input_image`: truecolor RGB image
- `landmask_binary`: binary landmask with 1=land, 0=water/ice

""" # TODO: add option to use alpha channel for mask
function apply_landmask(input_image::AbstractMatrix, landmask_binary::AbstractArray{Bool})
    image_masked = (.!landmask_binary) .* input_image
    return image_masked
end

function apply_landmask!(input_image::AbstractMatrix, landmask_binary::AbstractArray{Bool})
    input_image .= (.!landmask_binary) .* input_image
    return nothing
end
