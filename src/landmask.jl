"""
    create_landmask(landmask_image, struct_elem, fill_value_lower, fill_value_upper)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix, including a buffer to extend the land over any soft ice regions; land = 0, water/ice = 1. Retuns a named tuple with the dilated and non-dilated landmask.

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
    fill_value_upper::Int=2000
) where {T<:AbstractMatrix}
    landmask_binary = binarize_landmask(landmask_image)
    dilated = IceFloeTracker.MorphSE.dilate(landmask_binary, centered(struct_elem))
    return (dilated=ImageMorphology.imfill(.!dilated, (fill_value_lower, fill_value_upper)), non_dilated=landmask_binary)
end

function create_landmask(landmask_image)
    return create_landmask(landmask_image, make_landmask_se())
end

"""
    binarize_landmask(landmask_image)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix; land = 0, water/ice = 1.

# Arguments
- `landmask_image`: RGB land mask image from `fetchdata`
"""
function binarize_landmask(landmask_image::T)::BitMatrix where {T<:AbstractMatrix}
    if !(typeof(landmask_image) <: AbstractMatrix{Bool})
        landmask_no_dilate = Gray.(landmask_image) .> 0
    end
    return landmask_no_dilate
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

# in-place version
function apply_landmask!(input_image::AbstractMatrix, landmask_binary::BitMatrix)
    input_image .= landmask_binary .* input_image
    return nothing
end

"""
    remove_landmask(landmask, ice_mask)

Find the pixel indexes that are floating ice rather than soft or land ice. Returns an array of pixel indexes. 

# Arguments
- `landmask`: bitmatrix landmask for region of interest
- `ice_mask`: bitmatrix with ones equal to ice, zeros otherwise

"""
## NOTE(tjd): This function is called in `find_ice_labels.jl`
function remove_landmask(landmask::BitMatrix, ice_mask::BitMatrix)::Array{Int64}
    land = IceFloeTracker.apply_landmask(ice_mask, landmask)
    return [i for i ∈ 1:length(land) if land[i]]
end
