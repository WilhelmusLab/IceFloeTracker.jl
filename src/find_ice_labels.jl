"""
    find_reflectance_peaks(reflectance_channel, possible_ice_threshold;)

Find histogram peaks in single channels of a reflectance image and return the second greatest peak. If needed, edges can be returned as the first object from `build_histogram`. Similarly, peak values can be returned as the second object from `findmaxima`.

# Arguments
- `reflectance_channel`: either band 2 or band 1 of false-color reflectance image
- `possible_ice_threshold`: threshold value used to identify ice if not found on first or second pass

"""
function find_reflectance_peaks(
    reflectance_channel::Matrix{Float64}; possible_ice_threshold::Float64=Float64(75 / 255)
)::Int64
    reflectance_channel[reflectance_channel .< possible_ice_threshold] .= 0 #75 / 255
    _, counts = ImageContrastAdjustment.build_histogram(reflectance_channel)
    locs, _ = Peaks.findmaxima(counts)
    locs = sort(locs; rev=true)
    second_peak = locs[2]
    return second_peak
end

"""
    find_ice_labels(reflectance_image, landmask; band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Locate the pixels of obvious ice from false color reflectance image. Returns a binary mask with ice floes contrasted from background. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `reflectance_image`: corrected reflectance false color image - bands [7,2,1]
- `landmask`: bitmatrix landmask for region of interest
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)

"""
function find_ice_labels(
    reflectance_image::Matrix{RGB{Float64}},
    landmask::BitMatrix;
    band_7_threshold::Float64=Float64(5 / 255),
    band_2_threshold::Float64=Float64(230 / 255),
    band_1_threshold::Float64=Float64(240 / 255),
    band_7_threshold_relaxed::Float64=Float64(10 / 255),
    band_1_threshold_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)::Vector{Int64}

    ## Make ice masks
    cv = channelview(reflectance_image)
    mask_ice_band_7 = cv[1, :, :] .< band_7_threshold #5 / 255
    mask_ice_band_2 = cv[2, :, :] .> band_2_threshold #230 / 255
    mask_ice_band_1 = cv[3, :, :] .> band_1_threshold #240 / 255
    ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_1
    ice_labels = remove_landmask(landmask, ice)
    println("Done with masks")

    ## Find obvious ice floes
    if sum(abs.(ice_labels)) == 0
        mask_ice_band_7 = cv[1, :, :] .< band_7_threshold_relaxed #10 / 255
        mask_ice_band_1 = cv[3, :, :] .> band_1_threshold_relaxed #190 / 255
        ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_1
        ice_labels = remove_landmask(landmask, ice)
        if sum(abs.(ice_labels)) == 0
            ref_image_band_2 = cv[2, :, :]
            ref_image_band_1 = cv[3, :, :]
            band_2_peak = find_reflectance_peaks(ref_image_band_2)
            band_1_peak = find_reflectance_peaks(ref_image_band_1)
            mask_ice_band_2 = cv[2, :, :] .> band_2_peak / 255
            mask_ice_band_1 = cv[3, :, :] .> band_1_peak / 255
            ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_1
            ice_labels = remove_landmask(landmask, ice)
        end
    end
    println("Done with ice labels")
    return ice_labels
end