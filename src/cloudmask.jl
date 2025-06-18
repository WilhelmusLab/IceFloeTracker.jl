function convert_to_255_matrix(img)::Matrix{Int}
    img_clamped = clamp.(img, 0.0, 1.0)
    return round.(Int, img_clamped * 255)
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
    modis_721::AbstractArray{T};
    prelim_threshold::Real=(110 / 255),
    band_7_threshold::Real=(200 / 255),
    band_2_threshold::Real=(190 / 255),
    ratio_lower::Real=0.0,
    ratio_upper::Real=0.75,
)::BitMatrix where {T<:Union{AbstractRGB,TransparentRGB}}
    modis_band02 = green.(modis_721)
    modis_band07 = red.(modis_721)

    clouds_view = modis_band07 .> prelim_threshold

    # First find all the pixels that meet threshold logic in band 7 (channel 1) and band 2 (channel 2)
    # Masking clouds and discriminating cloud-ice
    mask_b7 = modis_band07 .< band_7_threshold
    mask_b2 = modis_band02 .> band_2_threshold
    mask_b7b2 = mask_b7 .&& mask_b2

    # Next find pixels that meet both thresholds and mask them from band 7 (channel 1) and band 2 (channel 2)
    b7_masked = mask_b7b2 .* modis_band07
    b2_masked = mask_b7b2 .* modis_band02

    cloud_ice_ratio = float.(b7_masked) ./ float.(b2_masked)
    mask_cloud_ice = @. (ratio_lower <= cloud_ice_ratio < ratio_upper)
    not_cloud = mask_cloud_ice .|| .!clouds_view
    cloud = .!not_cloud

    return cloud
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
    cloud::AbstractArray{Bool},
    false_color_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
)
    red_channel_with_non_clouds_zeroed = red.(false_color_image) .* cloud
    return colorview(Gray, red_channel_with_non_clouds_zeroed)
end
