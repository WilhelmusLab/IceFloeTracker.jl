"""
    create_landmask(landmask_image, struct_elem, fill_value_lower, fill_value_upper)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix, including a buffer to extend the land over any soft ice regions; land = 0, water/ice = 1.

# Arguments
- `landmask_image`: land mask image
- `struct_elem`: structuring element for dilation (optional)
- `fill_value_lower`: fill holes having at least these many pixels (optional)
- `fill_value_upper`: fill holes having at most these many pixels (optional)

"""
function create_landmask(
    landmask_image::T,
    struct_elem::AbstractMatrix{Bool};
    fill_value_lower::Int=0,
    fill_value_upper::Int=2000,
)::BitMatrix where T<:AbstractMatrix

    # binarize if not Boolean
    if !(typeof(landmask_image) <: AbstractMatrix{Bool})
        landmask_image = Gray.(landmask_image) .> 0
    end
    dilated = IceFloeTracker.MorphSE.dilate(landmask_image, struct_elem)
    return ImageMorphology.imfill(.!dilated, (fill_value_lower, fill_value_upper))
end

function create_landmask(landmask_image)
    create_landmask(landmask_image, make_landmask_se())
end

"""
    apply_landmask(input_image, landmask_binary)

Zero out pixels in land and soft ice regions on truecolor image, return RGB image with zero for all three channels on land/soft ice.


# Arguments
- `input_image`: truecolor RGB image
- `landmask_binary`: binary landmask with 1=land, 0=water/ice 

"""
function apply_landmask(input_image::AbstractMatrix, landmask_binary::BitMatrix)
    image_masked = landmask_binary .* input_image
    return image_masked
end
