using Images # I think this can be taken out, I have it here for the RGB types
abstract type AbstractCloudMaskAlgorithm end

struct LopezAcostaCloudMask <: AbstractCloudMaskAlgorithm
    prelim_threshold::Float64
    band_7_threshold::Float64
    band_2_threshold::Float64
    ratio_lower::Float64
    ratio_offset::Float64
    ratio_upper::Float64
    
    # enforce all are between 0 and 1 inclusive
    function LopezAcostaCloudMask(
        prelim_threshold, band_7_threshold, band_2_threshold,
        ratio_lower, ratio_offset, ratio_upper
    )
        0 ≤ prelim_threshold ≤ 1 || error("$prelim_threshold must be between 0 and 1")
        0 ≤ band_7_threshold ≤ 1 || error("$band_7_threshold must be between 0 and 1")
        0 ≤ band_2_threshold ≤ 1 || error("$band_2_threshold must be between 0 and 1")
        0 ≤ ratio_lower ≤ 1 || error("$ratio_lower must be between 0 and 1")
        0 ≤ ratio_upper ≤ 1 || error("$ratio_upper must be between 0 and 1")
        return new(
            prelim_threshold, band_7_threshold, band_2_threshold,
            ratio_lower, ratio_upper, ratio_offset
        )
    end
end

# set defaults to match LSW2019
# and enable named arguments
function LopezAcostaCloudMask(;
        prelim_threshold::Float64=110/255.,
        band_7_threshold::Float64=200/255.,
        band_2_threshold::Float64=190/255.,
        ratio_lower::Float64=0.0,
        ratio_offset::Float64=0.0,
        ratio_upper::Float64=0.75
)
    LopezAcostaCloudMask(prelim_threshold,
                         band_7_threshold,
                         band_2_threshold,
                         ratio_lower,
                         ratio_offset,
                         ratio_upper)
end


# use functor notation to define a function using the parameter struct
function (f::LopezAcostaCloudMask)(img::AbstractArray{<:Union{AbstractRGB,TransparentRGB}})
    mask_cloud_ice, clouds_view = _get_masks(
        img;
        prelim_threshold=f.prelim_threshold,
        band_7_threshold=f.band_7_threshold,
        band_2_threshold=f.band_2_threshold,
        ratio_lower=f.ratio_lower,
        ratio_offset=0.0,
        ratio_upper=f.ratio_upper,
    )
    return clouds_view .&& .!mask_cloud_ice
end

# define internal function to apply the LopezAcosta et al. mask
function _get_masks(
    false_color_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}};
    prelim_threshold::Float64,
    band_7_threshold::Float64,
    band_2_threshold::Float64,
    ratio_lower::Float64,
    ratio_offset::Float64,
    ratio_upper::Float64,
)

    # this method assumes the images are MODIS 7-2-1 false color images
    false_color_image_b7 = red.(false_color_image)
    false_color_image_b2 = green.(false_color_image)

    # clouds_view is "true" if pixel is likely cloud
    clouds_view = false_color_image_b7 .> prelim_threshold

    # next select potential ice pixels to unmask
    mask_b7 = false_color_image_b7 .< band_7_threshold
    mask_b2 = false_color_image_b2 .> band_2_threshold

    # finally unmask only the pixels where b2 and b7 are inside the "box"
    # and between the two ratios
    b2_masked = false_color_image_b2 .* (mask_b2 .&& mask_b7)
    b7_masked = false_color_image_b7 .* (mask_b2 .&& mask_b7)
    b7_masked_float, b2_masked_float = [Float64.(m) for m in [b7_masked, b2_masked]]
    b7_greater_than_adjusted_b2_lower = @. b7_masked_float >=
        (b2_masked_float * ratio_lower)
    b7_less_than_adjusted_b2_upper = @. b7_masked_float <
        (b2_masked_float * (ratio_upper - ratio_offset))
    mask_cloud_ice = b7_greater_than_adjusted_b2_lower .&& b7_less_than_adjusted_b2_upper

    # Returning the two masks for facilitating testing other related workflows such as conditional adaptive histogram equalization
    return mask_cloud_ice, clouds_view
end

# Deprecated
# this function is now only used in the test of the conditional adaptive histogram
# it was also used in an earlier version of the cloud mask algorithm
function convert_to_255_matrix(img)::Matrix{Int}
    img_clamped = clamp.(img, 0.0, 1.0)
    return round.(Int, img_clamped * 255)
end

# Potential upgrade: add method to create a sequence of cloud masks
"""
    create_cloudmask(img; f=AbstractCloudMaskAlgorithm)

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

# Arguments
- `false_color_image`: corrected reflectance false color image - bands [7,2,1]
- `prelim_threshold`: threshold value used to identify clouds in band 7, N0f8(RGB intensity/255)
- `band_7_threshold`: threshold value used to identify cloud-ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify cloud-ice in band 2, N0f8(RGB intensity/255)
- `ratio_lower`: threshold value used to set lower ratio of cloud-ice in bands 7 and 2
- `ratio_upper`: threshold value used to set upper ratio of cloud-ice in bands 7 and 2
- `ratio_offset`: offset value used to adjust the upper ratio of cloud-ice in bands 7 and 2

1. Lopez-Acosta, R., Schodlok, M. P., & Wilhelmus, M. M. (2019). Ice Floe Tracker: An algorithm to automatically retrieve Lagrangian trajectories via feature matching from moderate-resolution visual imagery. Remote Sensing of Environment, 234(111406), 1–15. (https://doi.org/10.1016/j.rse.2019.111406)[https://doi.org/10.1016/j.rse.2019.111406]
2. Watkins, D.M., Kim, M., Paniagua, C., Divoll, T., Holland, J.G., Hatcher, S., Hutchings, J.K., and Wilhelmus, M.M. (in prep). Calibration and validation of the Ice Floe Tracker algorithm. 
"""
function create_cloudmask(
    false_color_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    f::AbstractCloudMaskAlgorithm=LopezAcostaCloudMask(),
)
    return f(false_color_image)
end


"""
    apply_cloudmask(false_color_image, cloudmask)

Zero out pixels containing clouds where clouds and ice are not discernable. Arguments should be of the same size.

# Arguments
- `img`: RGB, RGBA, or Gray image to be masked
- `cloudmask`: binary cloudmask with clouds = 1, else = 0
- `modify_channel_1`: optional keyword argument for RGB images. If true, set the first channel to 0 in the returned image.
"""
function apply_cloudmask(
    img::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    cloudmask::AbstractArray{Bool};
    modify_channel_1::Bool=false,
)
    masked_image = deepcopy(img)
    masked_image[cloudmask] .= 0.0
    modify_channel_1 && begin
        # Specialty application where Band 7 is set to 0          
        image_view = channelview(masked_image)
        cloudmasked_view = StackedView(
            zeroarray, @view(image_view[2, :, :]), @view(image_view[3, :, :])
        )
        cloudmasked_image_rgb = colorview(RGB, cloudmasked_view)
        return cloudmasked_image_rgb
    end
    return masked_image
end

function apply_cloudmask(
    img::AbstractArray{<:Union{AbstractRGB,TransparentRGB,Gray}},
    cloudmask::AbstractArray{Bool},
)
    masked_image = deepcopy(img)
    masked_image[cloudmask] .= 0.0
    return masked_image
end

function apply_cloudmask!(
    img::AbstractArray{<:Union{AbstractRGB,TransparentRGB,Gray}},
    cloudmask::AbstractArray{Bool},
)
    img[cloudmask] .= 0.0
end



# dmw: in the future, we may want the option to use "missing".
# dmw: used only in ice-water-discrimination. could be generalized. compare to similar method in the conditional adaptive histogram
# is this not equivalent to something like apply_cloudmask(red.(img), cloudmask)? also note that this function still has the clouds=0 
# sense for the mask.
function create_clouds_channel(
    cloudmask::AbstractArray{Bool}, false_color_image::Matrix{RGB{Float64}}
)::Matrix{Gray{Float64}}
    return Gray.(@view(channelview(cloudmask .* false_color_image)[1, :, :]))
end