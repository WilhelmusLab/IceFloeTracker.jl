"""
Borrowing image algorithm and image filter types from ImageBinarization.jl
Filters are image algorithms whose input and output are both images.
Image binarization algorithms produce a binary image.
"""
abstract type AbstractImageAlgorithm end
abstract type AbstractImageFilter <: AbstractImageAlgorithm end

# dmw: could use AbstractImageBinarizationAlgorithm instead, but figured this is more specialized
abstract type AbstractCloudMaskAlgorithm <: AbstractImageFilter end

# dmw: want a default version that calls the LopezAcosta2019 cloud mask
create_cloudmask(img,
                 f::AbstractCloudMaskAlgorithm,
                 args...; kwargs...) = f(img, args...; kwargs...)

# dmw: the idea is to make a struct that
struct LopezAcosta2019 <: AbstractCloudMaskAlgorithm
    false_color_image::AbstractArray
    prelim_threshold::AbstractFloat
    band_7_threshold::AbstractFloat
    band_2_threshold::AbstractFloat
    ratio_lower::AbstractFloat
    ratio_upper::AbstractFloat
    return_all_layers::AbstractBool
    
    function LopezAcosta2019(false_color_image::AbstractArray,
                             prelim_threshold::AbstractFloat,
                             band_7_threshold::AbstractFloat,
                             band_2_threshold::AbstractFloat,
                             ratio_lower::AbstractFloat,
                             ratio_upper::AbstractFloat,
                             return_all_layers::AbstractBool)
        
        # Calibrated reflectance at 2.1-2.15 μm band
        false_color_image_b7 = false_color_image[1, :, :]

        # Calibrated reflectance at 0.84-0.88 μm band
        false_color_image_b2 = false_color_image[2, :, :]

        opaque_cloud, all_cloud = _get_cloud_masks_la2019(
            false_color_image_b7,
            false_color_image_b2,
            prelim_threshold=prelim_threshold,
            band_7_threshold=band_7_threshold,
            band_2_threshold=band_2_threshold,
            ratio_lower=ratio_lower,
            ratio_upper=ratio_upper,
        )

        return_all_layers && return opaque_cloud, all_cloud
        return opaque_cloud
    end
end




# dmw: potentially not needed; or to move to Util
function convert_to_255_matrix(img)::Matrix{Int}
    img_clamped = clamp.(img, 0.0, 1.0)
    return round.(Int, img_clamped * 255)
end





# dmw: changed name because we use multiple types of masks in the package
function _get_cloud_masks_la2019(
    band_2_reflectance::AbstractArray;
    band_7_reflectance::AbstractArray;
    prelim_threshold::AbstractFloat=Float64(110 / 255),
    band_7_threshold::AbstractFloat=Float64(200 / 255),
    band_2_threshold::AbstractFloat=Float64(190 / 255),
    ratio_lower::Number=Float64(0),
    ratio_upper::AbstractFloat=0.75
    )::Tuple{BitMatrix,BitMatrix}

    @assert 0 <= prelim_threshold <= 1
    @assert 0 <= band_7_threshold <= 1
    @assert 0 <= band_2_threshold <= 1
     
    # Estimate of total cloud cover
    all_cloud = band_7_reflectance .> prelim_threshold

    # First find all the pixels that meet threshold logic in band 7 (channel 1) and band 2 (channel 2)
    # These pixels are marked as potentially visible sea ice
    mask_b7 = band_7_reflectance .< band_7_threshold
    mask_b2 = band_2_reflectance .> band_2_threshold
    mask_b7b2 = mask_b7 .&& mask_b2

    # We check whether the ratio is small enough. Using leq instead of div prevents divide-by-zero instability
    ratio_test = band_7_reflectance .<= ratio_upper .* band_2_reflectance

    # Thin clouds are the pixels that are included in the "all clouds" but are unmasked by the threshold tests
    # Hence we can remove those to get the "opaque clouds" category.
    opaque_cloud = @. all_cloud && !(ratio_test && mask_b7b2)
    
    return opaque_cloud, all_cloud
end

"""
    create_cloudmask(false_color_image; prelim_threshold, band_7_threshold, band_2_threshold, ratio_lower, ratio_upper)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix; clouds = 1, else = 0. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406. The algorithm partitions the band 2 and band 7 intensities with a piecewise linear function, defined
by a series of thresholds. Thresholds can be specified as integers or scaled integers. The ratio offset shifts the sloped line in the algorithm; it is not present in the original 2019 version and was added to have a better fit to the observations in the 2025 algorithm update.

# Arguments
- `false_color_image`: corrected reflectance false color image - bands [7,2,1]
- `prelim_threshold`: threshold value used to identify clouds in band 7 (default 110 or 110/255)
- `band_7_threshold`: threshold value used to identify thin clouds over ice in band 7 (default 200 or 200/255)
- `band_2_threshold`: threshold value used to identify thin clouds over ice in band 2 (default 190 or 190/255)
- `ratio_offset`: horizontal shift of band 2 (default 0)
- `ratio_upper`: threshold value used to set upper limit for band 7 to band 2 reflectance ratio. Low values are marked cloud-free or thin cloud.
- `return_all_layers`: option to return the opaque cloud and all cloud separately. Default `false` returns only opaque cloud.

"""


function create_cloudmask(
    false_color_image::AbstractArray{<:Color3};
    prelim_threshold::Number=Float64(110 / 255), # Note: Using Number so that int or float can be provided
    band_7_threshold::Number=Float64(200 / 255),
    band_2_threshold::Number=Float64(190 / 255),
    ratio_lower::Number=Float(0),
    ratio_upper::AbstractFloat=0.75,
    return_all_layers=false
) # Note: Can return either a BitMatrix or a tuple of BitMatrices

    if 1 <= prelim_threshold <= 255
        prelim_threshold = prelim_threshold / 255.
        band_7_threshold = band_7_threshold / 255.
        band_2_threshold = band_2_threshold / 255.
    end
        
    opaque_cloud, all_cloud = _get_cloud_masks(
        false_color_image,
        prelim_threshold=prelim_threshold,
        band_7_threshold=band_7_threshold,
        band_2_threshold=band_2_threshold,
        ratio_offset=ratio_offset,
        ratio_upper=ratio_upper,
    )

    return_all_layers && return opaque_cloud, all_cloud
    return opaque_cloud
end

"""
    apply_cloudmask(false_color_image, cloudmask)

Zero out pixels containing clouds where clouds and ice are not discernable. Arguments should be of the same size.
If a color image is provided, then the first channel is replaced with 0's as well.

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