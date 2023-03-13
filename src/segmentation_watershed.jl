"""
    watershed_ice_floes(intermediate_segmentation_image;)
Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl to further isolate ice floes, returning a binary segmentation mask indicating potential sparse boundaries of ice floes.
# Arguments
-`intermediate_segmentation_image`: binary cloudmasked and landmasked intermediate file from segmentation B, either `SegB.not_ice_bit` or `SegB.ice_intersect`
"""
function watershed_ice_floes(intermediate_segmentation_image::BitMatrix)::BitMatrix
    features = Images.feature_transform(.!intermediate_segmentation_image)
    distances = 1 .- Images.distance_transform(features)
    seg_mask = ImageSegmentation.hmin_transform(distances, 2)
    seg_mask_bool = seg_mask .> 0
    markers = Images.label_components(seg_mask_bool)
    segment = ImageSegmentation.watershed(distances, markers)
    labels = ImageSegmentation.labels_map(segment)
    borders = Images.isboundary(labels)
    return borders
end

"""
    watershed_product(watershed_B_ice_intersect, watershed_B_not_ice;)
Intersects the outputs of watershed segmentation on intermediate files from segmentation B, indicating potential sparse boundaries of ice floes.
# Arguments
- `watershed_B_ice_intersect`: binary segmentation mask from `watershed_ice_floes`
- `watershed_B_not_ice`: binary segmentation mask from `watershed_ice_floes`
"""
function watershed_product(watershed_B_ice_intersect::BitMatrix, watershed_B_not_ice::BitMatrix;)::BitMatrix

    ## Intersect the two watershed files
    watershed_intersect = watershed_B_ice_intersect .* watershed_B_not_ice
    return watershed_intersect
end