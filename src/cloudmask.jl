"""
    create_cloudmask(reflectance_image; prelim_threshold, band_7_threshold, band_2_threshold, ratio_lower, ratio_upper)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix; clouds = 0, else = 1. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `reflectance_image`: corrected reflectance false color image - bands [7,2,1]
- `prelim_threshold`: threshold value used to identify clouds in band 7, N0f8(RGB intensity/255)
- `band_7_threshold`: threshold value used to identify cloud-ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify cloud-ice in band 2, N0f8(RGB intensity/255)
- `ratio_lower`: threshold value used to set lower ratio of cloud-ice in bands 7 and 2
- `ratio_upper`: threshold value used to set upper ratio of cloud-ice in bands 7 and 2

"""
function create_cloudmask(
    ref_image::Matrix{RGB{Float64}};
    prelim_threshold::Float64=Float64(110 / 255),
    band_7_threshold::Float64=Float64(200 / 255),
    band_2_threshold::Float64=Float64(190 / 255),
    ratio_lower::Float64=0.0,
    ratio_upper::Float64=0.75,
)::BitMatrix
    println("Setting thresholds")
    ref_view = channelview(ref_image)
    ref_image_b7 = ref_view[1, :, :]
    clouds_view = ref_image_b7 .> prelim_threshold
    mask_b7 = ref_image_b7 .< band_7_threshold
    mask_b2 = ref_view[2, :, :] .> band_2_threshold
    # First find all the pixels that meet threshold logic in band 7 (channel 1) and band 2 (channel 2)
    println("Masking clouds and discriminating cloud-ice")

    mask_b7b2 = mask_b7 .&& mask_b2
    # Next find pixels that meet both thresholds and mask them from band 7 (channel 1) and band 2 (channel 2)
    b7_masked = mask_b7b2 .* ref_image_b7
    b2_masked = mask_b7b2 .* ref_view[2, :, :]
    cloud_ice = Float64.(b7_masked) ./ Float64.(b2_masked)
    mask_cloud_ice = @. cloud_ice >= ratio_lower .&& cloud_ice < ratio_upper
    println("Creating final cloudmask")
    cloudmask = mask_cloud_ice .|| .!clouds_view
    return cloudmask
end

"""
    apply_cloudmask(reflectance_image, cloudmask)

Zero out pixels containing clouds where clouds and ice are not discernable. Arguments should be of the same size.

# Arguments
- `reference_image`: corrected reflectance false color image - bands [7,2,1] or grayscale
- `cloudmask`: binary cloudmask with clouds = 0, else = 1

"""
function apply_cloudmask(
    ref_image::Matrix{RGB{Float64}}, cloudmask::AbstractArray{Bool}
)::Matrix{RGB{Float64}}
    masked_image = cloudmask .* ref_image
    image_view = channelview(masked_image)
    cloudmasked_view = StackedView(zeroarray, image_view[2, :, :], image_view[3, :, :])
    cloudmasked_image_rgb = colorview(RGB, cloudmasked_view)
    return cloudmasked_image_rgb
end

function apply_cloudmask(
    ref_image::Matrix{Gray{Float64}}, cloudmask::AbstractArray{Bool}
)::Matrix{Gray{Float64}}
    return Gray.(cloudmask .* ref_image)
end

function create_clouds_channel(
    cloudmask::AbstractArray{Bool}, ref_image::Matrix{RGB{Float64}}
)::Matrix{Gray{Float64}}
    return Gray.(channelview(cloudmask .* ref_image)[1, :, :])
end
