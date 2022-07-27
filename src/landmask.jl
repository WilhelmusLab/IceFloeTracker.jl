"""
    create_landmask(landmask_image, num_pixels_dilate, fill_value)

Convert a 3-channel RGB land mask image to a 1-channel binary matrix, including a buffer to extend the land over any soft ice regions; land = 0, water/ice = 1.

# Arguments
- `landmask_image`: land mask image
- `struct_elem`: structuring element for dilation
- `fill_value`: number of pixels used to fill holes in land mask

"""
function create_landmask(
    landmask_image::Matrix{RGB{N0f8}}, struct_elem::Matrix{Bool}; fill_value::Int=2000
)::BitMatrix
    lm_binary = Gray.(landmask_image) .== 0
    radius = Int.(ceil.(size(struct_elem) ./ 2)[1]) #assumes symmetry
    pad_size = Fill(1, (radius, radius))
    lm_binary = IceFloeTracker.add_padding(lm_binary, pad_size)
    println("Dilation with strel")
    @time lm_binary_dilated = ImageProjectiveGeometry.imdilate(.!lm_binary, struct_elem)
    lm_binary_dilated = IceFloeTracker.remove_padding(lm_binary_dilated, pad_size)
    println("Closing any holes in mask")
    landmask_bool = (lm_binary_dilated .< 0.5)
    @time landmask_bool_filled = ImageMorphology.imfill(landmask_bool, (0, fill_value))
    return landmask_bool_filled
end

"""
    apply_landmask(input_image, landmask_binary)

Zero out pixels in land and soft ice regions on truecolor image, return RGB image with zero for all three channels on land/soft ice.


# Arguments
- `input_image`: truecolor RGB image
- `landmask_binary`: binary landmask with 1=land, 0=water/ice 

"""
function apply_landmask(input_image::Matrix, landmask_binary::BitMatrix)::Matrix
    image_masked = landmask_binary .* input_image
    return image_masked
end
