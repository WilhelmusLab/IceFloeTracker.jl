"""
    discriminate_ice_water(truecolor_image, landmask, struct_elem; lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Some text here

# Arguments
- `truecolor_image`: input image in truecolor
- `landmask`: bitmatrix landmask for region of interest
- `nbins`: number of bins during histogram build
# verify then add all the defaulted params

"""
function discriminate_ice_water(
    #truecolor_image::Matrix,
    reflectance_image::Matrix,
    normalized_image::Matrix;
    #landmask::BitMatrix,
    #cloudmask::BitMatrix;
    pad_size::Real=50,
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
    ref_image = IceFloeTracker.add_padding(
        load(reflectance_test_image_file)[test_region...],
        Pad(:replicate, (pad_size, pad_size)),
    )
    clouds_channel_masked = IceFloeTracker.apply_landmask(
        clouds_channel, landmask_bitmatrix
    )
    ref_image_masked = IceFloeTracker.apply_landmask(ref_image, landmask_bitmatrix)
    ref_image_1_masked = IceFloeTracker.apply_landmask(
        Gray.(ref_image_view[1, :, :]), landmask_bitmatrix
    )

    # first define all of the image variations

    image_cropped = IceFloeTracker.remove_padding(
        normalized_image, Pad((pad_size, pad_size), (pad_size, pad_size))
    ) # output during image normalization, landmasked
    image_clouds = IceFloeTracker.remove_padding(
        clouds_channel_masked, Pad((pad_size, pad_size), (pad_size, pad_size))
    ) # output during cloudmask apply, landmasked 
    image_cloudless = IceFloeTracker.remove_padding(
        ref_image_1_masked, Pad((pad_size, pad_size), (pad_size, pad_size))
    ) # channel 1 from source reflectance image, landmasked
    image_floes = IceFloeTracker.remove_padding(
        ref_image_masked, Pad((pad_size, pad_size), (pad_size, pad_size))
    ) # source reflectance, landmasked

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

    _, yyvals = imhist(image_clouds)
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

    return Z, Z2
end
