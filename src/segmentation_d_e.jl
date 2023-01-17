"""
    watershed_ice_floes(intermediate_segmentation_image;)

Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl or segmentation_c.jl to further isolate ice floes, returning a binary segmentation mask indicating potential sparse boundaries of ice floes.

# Arguments
-`intermediate_segmentation_image`: binary cloudmasked and landmasked intermediate file from segmentation B or C, typically either `segmentation_b_not_ice_mask` or `segmented_c`

"""
function watershed_ice_floes(intermediate_segmentation_image::BitMatrix)::BitMatrix
    features = Images.feature_transform(.!intermediate_segmentation_image)
    distances = 1 .- Images.distance_transform(features)
    seg_mask = ImageSegmentation.hmin_transform(distances, 2)
    seg_mask_bool = seg_mask .< 1
    markers = Images.label_components(seg_mask_bool)
    segment = ImageSegmentation.watershed(distances, markers)
    labels = ImageSegmentation.labels_map(segment)
    borders = Images.isboundary(labels)
    return borders
end

## function alias segmentation_D is watershed on the `segmented_c` from segmentation_C 
const segmentation_D = watershed_ice_floes
## function alias segmentation_E is watershed on the `not_ice_mask` from segmentation_B
const segmentation_E = watershed_ice_floes

"""
    segmentation_D_E(watershed_D, watershed_E;)

Intersects the outputs of watershed segmentation on intermediate files from segmentation B and C, indicating potential sparse boundaries of ice floes.

# Arguments
- `watershed_D`: binary segmentation mask from `segmentation_D`
- `watershed_E`: binary segmentation mask from `segmentation_E`

"""
function segmentation_D_E(watershed_D::BitMatrix, watershed_E::BitMatrix;)::BitMatrix

    ## Intersect the two watershed files
    watershed_intersect = watershed_D .* watershed_E
    return watershed_intersect
end
