using TiffImages
using LocalFilters
using Images

function create_landmask(landmask_image, num_pixels_dilate::Int, num_pixels_closing::Int)
    landmask_image = dropdims(landmask_image, dims = 3)
    landmask_binary = Gray.(landmask_image) .== 0
    landmask_binary = LocalFilters.dilate(.!landmask_binary, num_pixels_dilate)
    landmask_binary = LocalFilters.closing(landmask_binary, num_pixels_closing)
    return landmask_binary
end

function apply_landmask(input_image, landmask_binary::BitArray)
    image_masked = .!landmask_binary .* input_image
    return image_masked
end