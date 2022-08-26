"""
    segmentation_A(reflectance_image, ice_water_discriminated_image, landmask, cloudmask,band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

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
    ice_water_discriminated_image::Matrix,
    cloudmask::BitMatrix;
    num_clusters::Real=4,
    fuzziness::Float64=1.01,
    min_opening_area::Real=50,
    fill_range::Tuple=(0, 50),
)::BitMatrix
    num_clusters = num_clusters
    fuzziness = fuzziness
    input_image = RGB.(ice_water_discriminated_image)

    segmentation_result = fuzzy_cmeans(input_image, num_clusters, fuzziness)
    segmented_ice =
        Gray.(
            segmentation_result.centers[1] *
            reshape(segmentation_result.weights[:, 2], axes(input_image))
        ) * 1.5
    segmented_ice_cloudmasked = segmented_ice .* cloudmask
    segmented_ice_cloudmasked[segmented_ice_cloudmasked .> 0.5] .= 1
    segmented_ice_cloudmasked[segmented_ice_cloudmasked .< 0.5] .= 0
    segmented_ice_cloudmasked = convert(BitMatrix, segmented_ice_cloudmasked)

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
