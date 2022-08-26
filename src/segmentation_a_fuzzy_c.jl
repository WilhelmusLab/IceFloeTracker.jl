"""
    segmentation_A(ice_water_discriminated_image, cloudmask; num_clusters, fuzziness, min_opening_area, fill_range)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix with ice floes contrasted from background. Returns an image segmented and processed. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `ice_water_discrimination_image`: output image from `ice-water-discrimination.jl`
- `cloudmask`: bitmatrix cloudmask for region of interest
- `num_clusters`: the number of desired clusters/segmentation groups
- `fuzziness`: threshold to determine how much fuzziness to use during clustering
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
