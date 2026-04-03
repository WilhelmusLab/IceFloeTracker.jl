import ..Morphology: se_disk50
import Images: Gray
import OffsetArrays: centered
import Images: ImageMorphology, dilate

"""
    make_landmask_se()

Create a structuring element for dilating the landmask.
"""
make_landmask_se = se_disk50

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
    landmask_binary = binarize_landmask(landmask_image)
    dilated = dilate(landmask_binary, centered(struct_elem))
    return (
        dilated=.!ImageMorphology.imfill(.!dilated, (fill_value_lower, fill_value_upper)),
        non_dilated=landmask_binary,
    )
end

function create_landmask(landmask_image; strel=make_landmask_se())
    return create_landmask(landmask_image, strel)
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

# function binarize_landmask(landmask_image; tol=0.1)::BitMatrix
#     return binarize_landmask(Gray.(landmask_image); tol=tol)
# end

"""
    apply_mask(img, mask)

Zero out pixels in `img` where `mask` is `true`. Returns a new array.

# Arguments
- `img`: image array (e.g., RGB, Gray, BitMatrix, or any element type)
- `mask`: boolean mask; pixels where `mask` is `true` are set to zero
"""
function apply_mask(img::AbstractArray, mask::AbstractArray{Bool})
    masked_image = deepcopy(img)
    masked_image[mask] .= zero(eltype(img))
    return masked_image
end

"""
    apply_mask!(img, mask)

Zero out pixels in `img` where `mask` is `true`, modifying `img` in-place.

# Arguments
- `img`: image array (e.g., RGB, Gray, BitMatrix, or any element type)
- `mask`: boolean mask; pixels where `mask` is `true` are set to zero
"""
function apply_mask!(img::AbstractArray, mask::AbstractArray{Bool})
    img[mask] .= zero(eltype(img))
    return nothing
end
