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
- `possible_ice_threshold`: threshold value used to identify ice if not found on first or second pass

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
    if isempty(ice_labels)
        mask_ice_band_7 = cv[1, :, :] .< band_7_threshold_relaxed #10 / 255
        mask_ice_band_1 = cv[3, :, :] .> band_1_threshold_relaxed #190 / 255
        ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_1
        ice_labels = remove_landmask(landmask, ice)
        if isempty(ice_labels)
            ref_image_band_2 = cv[2, :, :]
            ref_image_band_1 = cv[3, :, :]
            ref_image_band_2[ref_image_band_2 .< possible_ice_threshold] .= 0 #75 / 255
            ref_image_band_1[ref_image_band_1 .< possible_ice_threshold] .= 0 #75 / 255
            edges_band_2, counts_band_2 = ImageContrastAdjustment.build_histogram(
                ref_image_band_2
            )
            locs_band_2, peaks_band_2 = Peaks.findmaxima(counts_band_2)
            locs_band_2 = sort(locs_band_2; rev=true)
            peak1 = locs_band_2[2]
            edges_band_1, counts_band_1 = ImageContrastAdjustment.build_histogram(
                ref_image_band_1
            )
            locs_band_1, peaks_band_1 = Peaks.findmaxima(counts_band_1)
            locs_band_1 = sort(locs_band_1; rev=true)
            peak2 = locs_band_1[2]
            mask_ice_band_2 = cv[2, :, :] .> peak1 / 255
            mask_ice_band_1 = cv[3, :, :] .> peak2 / 255
            ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_1
            ice_labels = remove_landmask(landmask, ice)
        end
    end
    println("Done with ice labels")
    return ice_labels
end
