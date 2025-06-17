function convert_to_255_matrix(img)::Matrix{Int}
    img_clamped = clamp.(img, 0.0, 1.0)
    return round.(Int, img_clamped * 255)
end

function _get_masks(
    false_color_image::Union{Matrix{RGB{Float64}}, Matrix{RGBA{N0f8}}};
    prelim_threshold::Float64=Float64(110 / 255),
    band_7_threshold::Float64=Float64(200 / 255),
    band_2_threshold::Float64=Float64(190 / 255),
    ratio_lower::Float64=0.0,
    ratio_upper::Float64=0.75,
    use_uint8::Bool=false,
)::Tuple{BitMatrix,BitMatrix}

    ref_view = channelview(false_color_image)
    false_color_image_b7 = @view ref_view[1, :, :]
    if use_uint8
        false_color_image_b7 = convert_to_255_matrix(false_color_image_b7)
    end

    clouds_view = false_color_image_b7 .> prelim_threshold
    mask_b7 = false_color_image_b7 .< band_7_threshold
    mask_b2 = @view(ref_view[2, :, :])
    if use_uint8
        mask_b2 = convert_to_255_matrix(mask_b2)
    end
    mask_b2 = mask_b2 .> band_2_threshold

    # First find all the pixels that meet threshold logic in band 7 (channel 1) and band 2 (channel 2)
    # Masking clouds and discriminating cloud-ice

    mask_b7b2 = mask_b7 .&& mask_b2

    # Next find pixels that meet both thresholds and mask them from band 7 (channel 1) and band 2 (channel 2)
    b7_masked = mask_b7b2 .* false_color_image_b7

    _b2 = @view(ref_view[2, :, :])
    b2_masked = use_uint8 ? convert_to_255_matrix(_b2) : _b2
    b2_masked = mask_b7b2 .* b2_masked

    cloud_ice = Float64.(b7_masked) ./ Float64.(b2_masked)
    mask_cloud_ice = @. (cloud_ice >= ratio_lower) && (cloud_ice < ratio_upper)

    return mask_cloud_ice, clouds_view
end


"""
    create_cloudmask(false_color_image; prelim_threshold, band_7_threshold, band_2_threshold, ratio_lower, ratio_upper)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix; clouds = 0, else = 1. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `false_color_image`: corrected reflectance false color image - bands [7,2,1]
- `prelim_threshold`: threshold value used to identify clouds in band 7, N0f8(RGB intensity/255)
- `band_7_threshold`: threshold value used to identify cloud-ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify cloud-ice in band 2, N0f8(RGB intensity/255)
- `ratio_lower`: threshold value used to set lower ratio of cloud-ice in bands 7 and 2
- `ratio_upper`: threshold value used to set upper ratio of cloud-ice in bands 7 and 2

"""
function create_cloudmask(
    false_color_image::Union{Matrix{RGB{Float64}}, Matrix{RGBA{N0f8}}};
    prelim_threshold::Float64=Float64(110 / 255),
    band_7_threshold::Float64=Float64(200 / 255),
    band_2_threshold::Float64=Float64(190 / 255),
    ratio_lower::Float64=0.0,
    ratio_upper::Float64=0.75,
)::BitMatrix
    mask_cloud_ice, clouds_view = _get_masks(
        false_color_image,
        prelim_threshold=prelim_threshold,
        band_7_threshold=band_7_threshold,
        band_2_threshold=band_2_threshold,
        ratio_lower=ratio_lower,
        ratio_upper=ratio_upper,
    )

    # Creating final cloudmask
    cloudmask = mask_cloud_ice .|| .!clouds_view
    return cloudmask
end

"""
    apply_cloudmask(false_color_image, cloudmask)

Zero out pixels containing clouds where clouds and ice are not discernable. Arguments should be of the same size.

# Arguments
- `false_color_image`: reference image, e.g. corrected reflectance false color image bands [7,2,1] or grayscale
- `cloudmask`: binary cloudmask with clouds = 0, else = 1

"""
function apply_cloudmask(
    false_color_image::Matrix{RGB{Float64}}, cloudmask::AbstractArray{Bool}
)::Matrix{RGB{Float64}}
    masked_image = cloudmask .* false_color_image
    image_view = channelview(masked_image)
    cloudmasked_view = StackedView(
        zeroarray, @view(image_view[2, :, :]), @view(image_view[3, :, :])
    )
    cloudmasked_image_rgb = colorview(RGB, cloudmasked_view)
    return cloudmasked_image_rgb
end

function apply_cloudmask(
    false_color_image::Matrix{Gray{Float64}}, cloudmask::AbstractArray{Bool}
)::Matrix{Gray{Float64}}
    return Gray.(cloudmask .* false_color_image)
end

function create_clouds_channel(
    cloudmask::AbstractArray{Bool}, false_color_image::Matrix{RGB{Float64}}
)::Matrix{Gray{Float64}}
    return Gray.(@view(channelview(cloudmask .* false_color_image)[1, :, :]))
end

