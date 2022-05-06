"""
    create_landmask(landmask_image, num_pixels_dilate, num_pixels_closing)

Convert a 3-channel RGB land mask TIFF to a 1-channel binary TIFF, including a buffer to extend the land over any soft ice regions; land = 0, water/ice = 1.

# Arguments
- `landmask_image`: loaded land mask TIFF from input
- `num_pixels_dilate`: number of pixels used to extend coast line; default = 50
- `num_pixels_closing`: number of pixels used to fill holes in land mask; default = 15

"""
function create_landmask(landmask_image::TiffImages.AbstractTIFF; num_pixels_dilate::Int=50, num_pixels_closing::Int=15)
    # Drop third dimension if it exists (test image had 3 dims: height x width x 1)
    landmask_image = dropdims(landmask_image, dims = 3)
    landmask_binary = Gray.(landmask_image) .== 0
    landmask_binary = LocalFilters.dilate(.!landmask_binary, num_pixels_dilate)
    landmask_binary = LocalFilters.closing(landmask_binary, num_pixels_closing)
    return landmask_binary
    # update to process inline
end

"""
    apply_landmask(input_image, landmask_binary)

Zero out pixels in land and soft ice regions on truecolor image, return RGB image with zero for all three channels on land/soft ice.


# Arguments
- `input_image`: truecolor RGB TIFF file stored in `input/truecolor`
- `landmask_binary`: binary landmask TIFF from `create_landmask`  

"""
function apply_landmask(input_image::TiffImages.AbstractTIFF, landmask_binary::BitArray)
    image_masked = .!landmask_binary .* input_image
    return image_masked
    # update to process inline
end