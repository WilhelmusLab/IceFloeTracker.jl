"""
    create_landmask(landmask_image, num_pixels_dilate, fill_value)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix, including a buffer to extend the land over any soft ice regions; land = 0, water/ice = 1.

# Arguments
- `landmask_image`: RGB or boolean image with land pixels represented as `true` or 1
- `struct_elem`: structuring element for dilation
- `fill_value`: number of pixels used to fill holes in land mask

"""
function create_landmask(
    landmask_image::Union{BitMatrix, Matrix{RGB{N0f8}}},
    struct_elem::Union{Matrix{Bool},BitMatrix};
    fill_value_lower::Int=0,
    fill_value_upper::Int=2000,
    )::BitMatrix

    # binarize image if not binarized
    if !(typeof(landmask_image) <: BitMatrix)
        lm_binary = BitArray(complement.(Gray.(landmask_image) .== 0)) # land pixels represented as ones. Faster than .! for inversion
    end
    
    println("Dilation with strel")
    # Dilate the land
    @time lm_binary_dilated = ImageMorphology.dilate(lm_binary, struct_elem);
    
    println("Closing any holes in mask")
    # ImageMorphology.imfill fills ones to zeros; input inversion required 
    @time landmask_bool_filled = ImageMorphology.imfill(
        .!(BitArray(lm_binary_dilated)), (fill_value_lower, fill_value_upper)) # faster than BitArray(complement.(lm_binary_dilated))
    
    # invert once more to get land pixels as ones for subsequent performant application
    return .!landmask_bool_filled; # faster than BitMatrix(complement.(img))
end

"""
    apply_landmask(input_image, landmask_binary)

Zero out pixels in land and soft ice regions on truecolor image, return RGB image with zero for all three channels on land/soft ice.


# Arguments
- `input_image`: truecolor RGB image
- `landmask_binary`: binary landmask with 1=land, 0=water/ice 

"""
function apply_landmask!(input_image::Matrix, landmask_binary::BitMatrix)
    return input_image[landmask_binary] .= 0
end
