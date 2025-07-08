function convert_to_255_matrix(img)::Matrix{Int}
    img_clamped = clamp.(img, 0.0, 1.0)
    return round.(Int, img_clamped * 255)
end

function _get_masks(
    false_color_image::Union{Matrix{RGB{Float64}},Matrix{RGBA{N0f8}},Matrix{RGB{N0f8}}};
    prelim_threshold::Float64=Float64(110 / 255),
    band_7_threshold::Float64=Float64(200 / 255),
    band_2_threshold::Float64=Float64(190 / 255),
    ratio_lower::Float64=0.0,
    ratio_upper::Float64=0.75,
    use_uint8::Bool=false,
    r_offset::Float64=0.0,
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

    b7_masked_float, b2_masked_float = [Float64.(m) for m in [b7_masked, b2_masked]]
    b7_greater_than_adjusted_b2_lower = @. b7_masked_float >=
        (b2_masked_float * ratio_lower)
    b7_less_than_adjusted_b2_upper = @. b7_masked_float <
        (b2_masked_float * (ratio_upper - r_offset))
    mask_cloud_ice = b7_greater_than_adjusted_b2_lower .&& b7_less_than_adjusted_b2_upper

    # Returning the two masks for facilitating testing other related workflows such as conditional adaptive histogram equalization
    return mask_cloud_ice, clouds_view
end

"""
    create_cloudmask(
    false_color_image::Union{Matrix{RGB{Float64}},Matrix{RGBA{N0f8}},Matrix{RGB{N0f8}}};
    prelim_threshold::Float64=Float64(110 / 255),
    band_7_threshold::Float64=Float64(200 / 255),
    band_2_threshold::Float64=Float64(190 / 255),
    ratio_lower::Float64=0.0,
    ratio_upper::Float64=0.75,
    r_offset::Float64=0.0,
)::BitMatrix

<<<<<<< HEAD
Convert a 3-channel false color reflectance image to a 1-channel binary matrix; clouds = 0, else = 1. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `false_color_image`: corrected reflectance false color image - bands [7,2,1]
- `prelim_threshold`: threshold value used to identify clouds in band 7, N0f8(RGB intensity/255)
- `band_7_threshold`: threshold value used to identify cloud-ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify cloud-ice in band 2, N0f8(RGB intensity/255)
- `ratio_lower`: threshold value used to set lower ratio of cloud-ice in bands 7 and 2
- `ratio_upper`: threshold value used to set upper ratio of cloud-ice in bands 7 and 2
- `r_offset`: offset value used to adjust the upper ratio of cloud-ice in bands 7 and 2

=======
Cloud masks in the IFT are BitMatrix objects such that for an image I and cloudmask C, cloudy pixels can be selected by I[C], and clear-sky pixels can be selected with I[.!C]. Construction of a cloud mask uses the syntax

```julia
f = CloudMaskAlgorithm(parameters)
C = create_cloudmask(img; CloudMaskAlgorithm)
```

By default, `create_cloudmask` uses the algorithm found in [1]. This algorithm converts a 3-channel MODIS 7-2-1 false color image into a 1-channel binary matrix in which clouds = 1 and anything else = 0. The algorithm aims to identify patches of opaque cloud while allowing thin and transparent cloud to remain. This algorithm is instantiated using

```julia
f = LopezAcostaCloudMask()
```

In this case, the default values are applied. It can also called using a set of customized parameters. These values must be real numbers between 0 and 1. To reproduce the default parameters, you may call

```julia
f = LopezAcostaCloudMask(prelim_threshold=110/255, band_7_threshold=200/255, band_2_threshold=190/255, ratio_lower=0.0, ratio_upper=0.75).
```

A stricter cloud mask was defined in [2], covering more cloudy pixels while minimally impacting the masking of cloud-covered ice pixels.

```julia
f = LopezAcostaCloudMask(prelim_threshold=53/255, band_7_threshold=130/255, band_2_threshold=169/255, ratio_lower=0.0, ratio_upper=0.53).
```

These parameters together define a piecewise linear partition of pixels based on their Band 7 and Band 2 callibrated reflectance. Pixels with intensity above `prelim_threshold` are considered as potential cloudy pixels. Then, pixels with Band 7 reflectance less than `band_7_threshold`, Band 2 reflectance greater than `band_2_threshold`, and Band 7 to Band 2 ratios between `ratio_lower` and `ratio_upper` are removed from the cloud mask (i.e., set to cloud-free).


1. Lopez-Acosta, R., Schodlok, M. P., & Wilhelmus, M. M. (2019). Ice Floe Tracker: An algorithm to automatically retrieve Lagrangian trajectories via feature matching from moderate-resolution visual imagery. Remote Sensing of Environment, 234(111406), 1–15. (https://doi.org/10.1016/j.rse.2019.111406)[https://doi.org/10.1016/j.rse.2019.111406]
2. Watkins, D.M., Kim, M., Paniagua, C., Divoll, T., Holland, J.G., Hatcher, S., Hutchings, J.K., and Wilhelmus, M.M. (in prep). Calibration and validation of the Ice Floe Tracker algorithm. 
>>>>>>> c34e38f (formatting docstring)
"""
function create_cloudmask(
    false_color_image::Union{Matrix{RGB{Float64}},Matrix{RGBA{N0f8}},Matrix{RGB{N0f8}}};
    prelim_threshold::Float64=Float64(110 / 255),
    band_7_threshold::Float64=Float64(200 / 255),
    band_2_threshold::Float64=Float64(190 / 255),
    ratio_lower::Float64=0.0,
    ratio_upper::Float64=0.75,
    r_offset::Float64=0.0,
)::BitMatrix
    mask_cloud_ice, clouds_view = _get_masks(
        false_color_image;
        prelim_threshold=prelim_threshold,
        band_7_threshold=band_7_threshold,
        band_2_threshold=band_2_threshold,
        ratio_lower=ratio_lower,
        ratio_upper=ratio_upper,
        r_offset=r_offset,
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
    false_color_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    cloudmask::AbstractArray{Bool},
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
