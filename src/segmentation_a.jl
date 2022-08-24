"""
    segmentation_A(reflectance_image, ice_water_discriminated_image, landmask, band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix with ice floes contrasted from background. Returns an image segmented and processed. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `ref_image`: corrected reflectance false color image - bands [7,2,1]
- `ice_water_discrimination_image`: output image from `ice-water-discrimination.jl`
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)
- `possible_ice_threshold`: threshold value used to identify ice if not found on first or second pass

"""
function remove_landmask(landmask::BitMatrix, ice_mask::BitMatrix)::Array
    data_no_landmask = []
    data = IceFloeTracker.apply_landmask(convert(Matrix, ice_mask), landmask)
    for (idx, val) in enumerate(data)
        if val != 0
            push!(data_no_landmask, idx)
        end
    end
    return data_no_landmask
end

function segmentation_A()::BitMatrix
    ice_water_discriminated_image = Array{Float32,2}(ice_water_discriminated_image)
    height, width = size(ice_water_discriminated_image)
    data = reshape(ice_water_discriminated_image, 1, height * width)
    classes = kmeans(data, 4; maxiter=50, display=:iter, init=:kmpp)
    class_assignments = assignments(classes)
    segmented = Gray.(((reshape(class_assignments, height, width)) .- 1) ./ 3) ## pixel_labels in matlab

    ## Make ice masks
    landmask_bitmatrix = convert(BitMatrix, load(current_landmask_file))
    ref_image = load(reflectance_test_image_file)[test_region...]
    cv = channelview(ref_image)
    mask_ice_band_7 = cv[1, :, :] .< band_7_threshold #5 / 255
    mask_ice_band_2 = cv[2, :, :] .> band_2_threshold #230 / 255
    mask_ice_band_1 = cv[3, :, :] .> band_1_threshold #240 / 255
    ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_3
    ice_labels = remove_landmask(landmask_bitmatrix, ice)

    if isempty(ice_labels)
        mask_ice_band_7 = cv[1, :, :] .< band_7_threshold_relaxed #10 / 255
        mask_ice_band_1 = cv[3, :, :] .> band_1_threshold_relaxed #190 / 255
        ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_3
        ice_labels = remove_landmask(landmask_bitmatrix, ice)
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
            ice = mask_ice_1 .&& mask_ice_2 .&& mask_ice_3
            ice_labels = remove_landmask(landmask_bitmatrix, ice)
            nlabel = mode(segmented[ice_labels])
        else
            nlabel = mode(segmented[ice_labels])
        end
    else
        nlabel = mode(segmented[ice_labels])
    end

    segmented_ice = segmented .== nlabel

    segmented_ice_cloudmasked = segmented_ice .* cloudmask

    segmented_ice_opened = ImageMorphology.area_opening(
        segmented_ice_cloudmasked; min_area=80
    ) #try diff values here; BW_test

    segmented_ice_filled =
        Gray.(ImageMorphology.imfill(convert(BitMatrix, (segmented_ice_opened)), (0, 200))) #BW_test3

    segmented_ice_filled_comp = Gray.(.!segmented_ice_filled)

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled_comp

    segmented_A = Gray.(segmented_ice_cloudmasked .|| diff_matrix)

    return segmented_A
end
