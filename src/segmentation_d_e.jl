"""
    watershed_ice_floes(intermediate_segmentation_image;)

Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl or segmentation_c.jl to further isolate ice floes, returning a binary segmentation mask indicating potential sparse boundaries of ice floes.

# Arguments
-`intermediate_segmentation_image`: binary cloudmasked and landmasked intermediate file from segmentation B or C, typically either `segmentation_b_not_ice_mask` or `segmented_c`

"""
function watershed_ice_floes(intermediate_segmentation_image::BitMatrix)::BitMatrix
    features = ImageSegmentation.feature_transform(.!intermediate_segmentation_image)
    distances = -1 .* ImageSegmentation.distance_transform(features)
    seg_mask = ImageSegmentation.hmin_transform(distances, 2)
    seg_mask_minima = ImageSegmentation.local_minima(seg_mask; connectivity=2)
    seg_mask_minima[seg_mask_minima .> 0] .= 1
    seg_mask_bool = Bool.(seg_mask_minima)
    markers = ImageSegmentation.label_components(seg_mask_minima)
    segment = ImageSegmentation.watershed(distances, markers; mask=seg_mask_bool)
    labels = ImageSegmentation.labels_map(segment)
    watershed_bitmatrix = labels .!= 0

    return watershed_bitmatrix
end

## function alias segmentation_D is watershed on the `not_ice_mask` from segmentation_B
const segmentation_D = watershed_ice_floes
## function alias segmentation_E is watershed on the `segmented_c` from segmentation_C
const segmentation_E = watershed_ice_floes

"""
    segmentation_D_E(watershed_B, watershed_C;)

Intersects the outputs of watershed segmentation on intermediate files from segmentation B and C, indicating potential sparse boundaries of ice floes.

# Arguments
- `watershed_B`: binary cloudmasked and landmasked segmentation mask from `segmentation_D`
- `watershed_C`: binary cloudmasked and landmasked segmentation mask from `segmentation_E`

"""
function segmentation_D_E(watershed_B::BitMatrix, watershed_C::BitMatrix;)::BitMatrix

    ## Intersect the two watershed files
    watershed_intersect = watershed_B .* watershed_C

    return watershed_intersect
end
