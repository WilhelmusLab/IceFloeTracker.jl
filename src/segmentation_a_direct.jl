"""
    remove_landmask(landmask, ice_mask)

Find the pixel indexes that are floating ice rather than soft or land ice. Returns an array of pixel indexes. 

# Arguments
- `landmask`: bitmatrix landmask for region of interest
- `ice_mask`: bitmatrix with ones equal to ice, zeros otherwise

"""
function remove_landmask(landmask::BitMatrix, ice_mask::BitMatrix)::Array{Int64}
    indexes_no_landmask = []
    land = IceFloeTracker.apply_landmask(ice_mask, landmask)
    for (idx, val) in enumerate(land)
        if val != 0
            push!(indexes_no_landmask, idx)
        end
    end
    return indexes_no_landmask
end

"""
    segmentation_A(reflectance_image, ice_water_discriminated_image, landmask, cloudmask;band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold, min_opening_area, fill_range)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix with ice floes contrasted from background. Returns an image segmented and processed. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `reflectance_image`: corrected reflectance false color image - bands [7,2,1]
- `ice_water_discrimination_image`: output image from `ice-water-discrimination.jl`
- `landmask`: bitmatrix landmask for region of interest
- `cloudmask`: bitmatrix cloudmask for region of interest
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)
- `possible_ice_threshold`: threshold value used to identify ice if not found on first or second pass
- `min_opening_area`: minimum size of pixels to use during morphoilogical opening
- `fill_range`: range of values dictating the size of holes to fill

"""
function segmentation_A(
    reflectance_image::Matrix{RGB{Float64}},
    ice_water_discriminated_image::Matrix{Gray{Float64}},
    landmask::BitMatrix,
    cloudmask::BitMatrix;
    band_7_threshold::Float64=Float64(5 / 255),
    band_2_threshold::Float64=Float64(230 / 255),
    band_1_threshold::Float64=Float64(240 / 255),
    band_7_threshold_relaxed::Float64=Float64(10 / 255),
    band_1_threshold_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
    min_opening_area::Real=50,
    fill_range::Tuple=(0, 50),
)::BitMatrix
    ice_water_discriminated_image = Matrix{Float64}(ice_water_discriminated_image)
    ice_water_discrimination_height, ice_water_discrimination_width = size(
        ice_water_discriminated_image
    )
    ice_water_discriminateed_1d = reshape(
        ice_water_discriminated_image,
        1,
        ice_water_discrimination_height * ice_water_discrimination_width,
    )
    feature_classes = Clustering.kmeans(
        ice_water_discriminateed_1d, 4; maxiter=50, display=:iter, init=:kmpp
    )
    class_assignments = assignments(feature_classes)

    ## NOTE(tjd): this reshapes column major vector of kmeans classes back into original image shape
    segmented = reshape(
        class_assignments, ice_water_discrimination_height, ice_water_discrimination_width
    )

    ## Make ice masks
    cv = channelview(reflectance_image)
    mask_ice_band_7 = cv[1, :, :] .< band_7_threshold #5 / 255
    mask_ice_band_2 = cv[2, :, :] .> band_2_threshold #230 / 255
    mask_ice_band_1 = cv[3, :, :] .> band_1_threshold #240 / 255
    ice = mask_ice_band_7 .&& mask_ice_band_2 .&& mask_ice_band_1
    ice_labels = remove_landmask(landmask, ice)

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
            nlabel = StatsBase.mode(segmented[ice_labels])
        else
            nlabel = StatsBase.mode(segmented[ice_labels])
        end
    else
        nlabel = StatsBase.mode(segmented[ice_labels])
    end

    ## Isolate ice floes and contrast from background
    segmented_ice = segmented .== nlabel
    segmented_ice_cloudmasked = segmented_ice .* cloudmask

    segmented_ice_opened = ImageMorphology.area_opening(
        segmented_ice_cloudmasked; min_area=min_opening_area
    ) #BW_test in matlab code

    segmented_opened_flipped = .!segmented_ice_opened

    segmented_ice_filled = ImageMorphology.imfill(
        convert(BitMatrix, (segmented_opened_flipped)), fill_range
    ) #BW_test3 in matlab code

    segmented_ice_filled_comp = complement.(segmented_ice_filled)

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled_comp

    segmented_A = segmented_ice_cloudmasked .|| diff_matrix

    return segmented_A
end
