"""
    discriminate_ice_water(reflectance_image, reflectance_image_band7, normalized_image, landmask_bitmatrix, clouds_channel; floes_threshold, mask_clouds_lower, mask_clouds_upper, kurt_thresh_lower, kurt_thresh_upper, skew_thresh, st_dev_thresh_lower, st_dev_thresh_upper, clouds2_threshold, differ_threshold, nbins)

Generates an image with ice floes apparent after filtering and combining previously processed versions of reflectance and truecolor images from the same region of interest. Returns an image ready for segmentation to isolate floes.

# Arguments
- `reflectance_image`: input image in false color reflectance
- `reflectance_image_band7`: first channel (band 7) of reflectance image, output from `cloudmask.jl`
- `normalized_image`: a normalized, landmasked truecolor image, output from `normalization.jl`
- `landmask_bitmatrix`: bitmatrix landmask for region of interest
- `clouds_channel`: first channel of cloudmasked reflectance image, output from `cloudmask.jl`
- `floes_threshold`: heuristic applied to original reflectance image
- `mask_clouds_lower`: lower heuristic applied to mask out clouds
- `mask_clouds_upper`: upper heuristic applied to mask out clouds
- kurt_thresh_lower`: lower heuristic used to set pixel value threshold based on kurtosis in histogram
- `kurt_thresh_upper`: upper heuristic used to set pixel value threshold based on kurtosis in histogram
- `skew_thresh`: heuristic used to set pixel value threshold based on skewness in histogram
- `st_dev_thresh_lower`: lower heuristic used to set pixel value threshold based on standard deviation in histogram
- `st_dev_thresh_upper`: upper heuristic used to set pixel value threshold based on standard deviation in histogram
- `clouds2_threshold`: heuristic used to set pixel value threshold based on ratio of clouds
- `differ_threshold`: heuristic used to calculate proportional intensity in histogram
- `nbins`: number of bins during histogram build
# verify then add all the defaulted params

"""
function discriminate_ice_water(
    reflectance_image::Matrix,
    reflectance_image_band7::Matrix,
    normalized_image::Matrix,
    landmask_bitmatrix::BitMatrix,
    clouds_channel::Matrix;
    floes_threshold::N0f8=N0f8(100 / 255),
    mask_clouds_lower::N0f8=N0f8(17 / 255),
    mask_clouds_upper::N0f8=N0f8(30 / 255),
    kurt_thresh_lower::Real=2,
    kurt_thresh_upper::Real=8,
    skew_thresh::Real=4,
    st_dev_thresh_lower::N0f8=N0f8(84 / 255),
    st_dev_thresh_upper::N0f8=N0f8(98.9 / 255),
    clouds2_threshold::Float64=0.02,
    differ_threshold::Float64=0.6,
    nbins::Real=155,
)::Matrix
    # first define all of the image variations

    image_cropped = normalized_image # output during image normalization, landmasked
    image_clouds = IceFloeTracker.apply_landmask(clouds_channel, landmask_bitmatrix) # output during cloudmask apply, landmasked 
    image_cloudless = IceFloeTracker.apply_landmask(
        reflectance_image_band7, landmask_bitmatrix
    ) # channel 1 (band 7) from source reflectance image, landmasked
    image_floes = IceFloeTracker.apply_landmask(reflectance_image, landmask_bitmatrix) # source reflectance, landmasked

    image_floes_view = channelview(image_floes)

    floes_band_2 = image_floes_view[2, :, :]
    floes_band_1 = image_floes_view[3, :, :]

    floes_band_2_keep = floes_band_2[floes_band_2 .> floes_threshold]

    floes_band_1_keep = floes_band_1[floes_band_1 .> floes_threshold]

    _, yyy = ImageContrastAdjustment.build_histogram(floes_band_2_keep, nbins)

    _, vals = Peaks.findmaxima(yyy)

    if isempty(vals)
        Z2 = zeroarray(size(image_floes))
    else
        differ = vals / (maximum(vals))
        proportional_intensity = sum(differ .> differ_threshold) / length(differ)
    end

    kurt_band_2 = kurtosis(floes_band_2_keep)
    skew_band_2 = skewness(floes_band_2_keep)
    kurt_band_1 = kurtosis(floes_band_1_keep)

    standard_dev = std(Float64.(image_cropped))

    _, yyvals = build_histogram(image_clouds)
    clouds1 = sum(yyvals[51:end])
    total1 = sum(yyvals)
    clouds2 = clouds1 / total1

    threshold_50_check = (
        (
            (abs(kurt_band_2 > kurt_thresh_upper)) ||
            (abs(kurt_band_2 < kurt_thresh_lower)) &&
            ((abs(kurt_band_1 > kurt_thresh_upper)))
        ) ||
        (
            (abs(kurt_band_2 < kurt_thresh_lower)) &&
            (abs(skew_band_2 < skew_thresh)) &&
            proportional_intensity < 0.1
        ) ||
        proportional_intensity < 0.01
    )
    threshold_130_check =
        (clouds2 .< clouds2_threshold && standard_dev > st_dev_thresh_lower) ||
        (standard_dev > st_dev_thresh_upper)

    if threshold_50_check
        THRESH = 50 / 255
    elseif threshold_130_check
        THRESH = 130 / 255
    else
        THRESH = 80 / 255
    end

    D1 = copy(image_cropped)
    D1[D1 .> THRESH] .= 0
    Z = image_cropped - (D1 * 3)

    D2 = copy(image_cloudless)
    mask_image_clouds = (
        image_clouds .< mask_clouds_lower .|| image_clouds .> mask_clouds_upper
    )
    D2 = D2 .* .!mask_image_clouds
    Z2 = Z - (D2 * 3)

    return Z2
end
