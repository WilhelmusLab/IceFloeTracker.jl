"""
    discriminate_ice_water(
    falsecolor_image::Matrix{RGB{Float64}},
    normalized_image::Matrix{Gray{Float64}},
    landmask_bitmatrix::T,
    cloudmask_bitmatrix::T,
    floes_threshold::Float64=Float64(100 / 255),
    mask_clouds_lower::Float64=Float64(17 / 255),
    mask_clouds_upper::Float64=Float64(30 / 255),
    kurt_thresh_lower::Real=2,
    kurt_thresh_upper::Real=8,
    skew_thresh::Real=4,
    st_dev_thresh_lower::Float64=Float64(84 / 255),
    st_dev_thresh_upper::Float64=Float64(98.9 / 255),
    clouds_ratio_threshold::Float64=0.02,
    differ_threshold::Float64=0.6,
    nbins::Real=155
)



Generates an image with ice floes apparent after filtering and combining previously processed versions of falsecolor and truecolor images from the same region of interest. Returns an image ready for segmentation to isolate floes.

Note: This function mutates the landmask object to avoid unnecessary memory allocation. If you need the original landmask, make a copy before passing it to this function. Example: `discriminate_ice_water(falsecolor_image, normalized_image, copy(landmask_bitmatrix), cloudmask_bitmatrix)`

# Arguments
- `falsecolor_image`: input image in false color reflectance
- `normalized_image`: normalized version of true color image
- `landmask_bitmatrix`: landmask for region of interest
- `cloudmask_bitmatrix`: cloudmask for region of interest
- `floes_threshold`: heuristic applied to original false color image
- `mask_clouds_lower`: lower heuristic applied to mask out clouds
- `mask_clouds_upper`: upper heuristic applied to mask out clouds
- `kurt_thresh_lower`: lower heuristic used to set pixel value threshold based on kurtosis in histogram
- `kurt_thresh_upper`: upper heuristic used to set pixel value threshold based on kurtosis in histogram
- `skew_thresh`: heuristic used to set pixel value threshold based on skewness in histogram
- `st_dev_thresh_lower`: lower heuristic used to set pixel value threshold based on standard deviation in histogram
- `st_dev_thresh_upper`: upper heuristic used to set pixel value threshold based on standard deviation in histogram
- `clouds2_threshold`: heuristic used to set pixel value threshold based on ratio of clouds
- `differ_threshold`: heuristic used to calculate proportional intensity in histogram
- `nbins`: number of bins during histogram build

"""
function discriminate_ice_water(
    falsecolor_image::AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    normalized_image::AbstractMatrix{<:Gray},
    landmask_bitmatrix::T,
    cloudmask_bitmatrix::T;
    floes_threshold::Real=(100 / 255),
    mask_clouds_lower::Real=(17 / 255),
    mask_clouds_upper::Real=(30 / 255),
    kurt_thresh_lower::Real=2,
    kurt_thresh_upper::Real=8,
    skew_thresh::Real=4,
    st_dev_thresh_lower::Real=(84 / 255),
    st_dev_thresh_upper::Real=(98.9 / 255),
    clouds_ratio_threshold::Real=0.02,
    differ_threshold::Real=0.6,
    nbins::Real=155,
)::AbstractMatrix where {T<:AbstractArray{Bool}}
    clouds_channel = IceFloeTracker.create_clouds_channel(
        cloudmask_bitmatrix, falsecolor_image
    )
    falsecolor_image_band7 = red.(falsecolor_image)

    # first define all of the image variations
    image_clouds = IceFloeTracker.apply_landmask(clouds_channel, landmask_bitmatrix) # output during cloudmask apply, landmasked
    image_cloudless = IceFloeTracker.apply_landmask(
        falsecolor_image_band7, landmask_bitmatrix
    ) # channel 1 (band 7) from source falsecolor image, landmasked
    image_floes = IceFloeTracker.apply_landmask(falsecolor_image, landmask_bitmatrix) # source false color reflectance, landmasked

    floes_band_2 = green.(image_floes)
    floes_band_1 = blue.(image_floes)

    # keep pixels greater than intensity 100 in bands 2 and 1
    floes_band_2_keep = floes_band_2[floes_band_2 .> floes_threshold]
    floes_band_1_keep = floes_band_1[floes_band_1 .> floes_threshold]

    _, floes_bin_counts = ImageContrastAdjustment.build_histogram(floes_band_2_keep, nbins)
    _, vals = Peaks.findmaxima(floes_bin_counts)

    differ = vals / (maximum(vals))
    proportional_intensity = sum(differ .> differ_threshold) / length(differ) # finds the proportional intensity of the peaks in the histogram

    # compute kurtosis, skewness, and standard deviation to use in threshold filtering
    kurt_band_2 = kurtosis(floes_band_2_keep)
    skew_band_2 = skewness(floes_band_2_keep)
    kurt_band_1 = kurtosis(floes_band_1_keep)
    standard_dev = stdmult(⋅, normalized_image)

    # find the ratio of clouds in the image to use in threshold filtering
    _, clouds_bin_counts = ImageContrastAdjustment.build_histogram(image_clouds .> 0)
    total_clouds = sum(clouds_bin_counts[51:end])
    total_all = sum(clouds_bin_counts)
    clouds_ratio = total_clouds / total_all

    threshold_50_check = _check_threshold_50(
        kurt_band_1,
        kurt_band_2,
        kurt_thresh_lower,
        kurt_thresh_upper,
        skew_band_2,
        skew_thresh,
        proportional_intensity,
    )

    threshold_130_check = _check_threshold_130(
        clouds_ratio,
        clouds_ratio_threshold,
        standard_dev,
        st_dev_thresh_lower,
        st_dev_thresh_upper,
    )

    if threshold_50_check
        THRESH = 50 / 255
    elseif threshold_130_check
        THRESH = 130 / 255
    else
        THRESH = 80 / 255 #intensity value of 80
    end

    normalized_image_copy = copy(normalized_image)
    normalized_image_copy[normalized_image_copy .> THRESH] .= 0
    @. normalized_image_copy = normalized_image - (normalized_image_copy * 3)

    mask_image_clouds = (
        image_clouds .< mask_clouds_lower .|| image_clouds .> mask_clouds_upper
    )

    band7_masked = image_cloudless .* mask_image_clouds

    ice_water_discriminated_image =
        clamp01nan.(normalized_image_copy .- (band7_masked .* 3))

    return ice_water_discriminated_image
end

function _check_threshold_50(
    kurt_band_1,
    kurt_band_2,
    kurt_thresh_lower,
    kurt_thresh_upper,
    skew_band_2,
    skew_thresh,
    proportional_intensity,
)
    return ( # intensity value of 50
        (
            (kurt_band_2 > kurt_thresh_upper) ||
            (kurt_band_2 < kurt_thresh_lower) && (kurt_band_1 > kurt_thresh_upper)
        ) ||
        (
            (kurt_band_2 < kurt_thresh_lower) &&
            (skew_band_2 < skew_thresh) &&
            proportional_intensity < 0.1
        ) ||
        proportional_intensity < 0.01
    )
end

function _check_threshold_130(
    clouds_ratio,
    clouds_ratio_threshold,
    standard_dev,
    st_dev_thresh_lower,
    st_dev_thresh_upper,
)
    return (clouds_ratio .< clouds_ratio_threshold && standard_dev > st_dev_thresh_lower) ||
           (standard_dev > st_dev_thresh_upper)
end
