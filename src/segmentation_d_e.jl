"""
    segmentation_D(segmentation_b_not_ice_mask, segmented_c;)

Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl and segmentation_c.jl to further isolate ice floes, returning two binary segmentation masks.

# Arguments
- `segmentation_b_not_ice_mask`: binary cloudmasked and landmasked intermediate file `not_ice_mask` from `segmentation_b.jl`
- `segmented_c`: binary cloudmasked and landmasked output file from `segmentation_c.jl`


"""
function segmentation_D_E(
    segmentation_b_not_ice_mask::BitMatrix, segmented_c::BitMatrix;
)::Tuple{BitMatrix,BitMatrix,BitMatrix}

    ## Watershed on segmentation_B intermediate
    features_B = feature_transform(.!segmentation_b_not_ice_mask)
    distances_B = 1 .- (distance_transform(features_B))
    seg_B_mask = ImageSegmentation.hmin_transform(distances_B, 2)
    seg_B_mask_minima = local_minima(seg_B_mask)
    seg_B_mask_minima[seg_B_mask_minima .> 0] .= 1
    seg_B_mask_bool = Bool.(seg_B_mask_minima)
    markers_B = label_components(seg_B_mask_minima)
    segment_B = watershed(distances_B, markers_B; mask=seg_B_mask_bool)
    labels_B = labels_map(segment_B)
    watershed_B = ifelse.(labels_B .== 0, 0, 1)

    ## Watershed on segmentation_C
    features_C = feature_transform(.!segmented_c)
    distances_C = 1 .- (distance_transform(features_C))
    seg_C_mask = ImageSegmentation.hmin_transform(distances_C, 2)
    seg_C_mask_minima = local_minima(seg_C_mask)
    seg_C_mask_minima[seg_C_mask_minima .> 0] .= 1
    seg_C_mask_bool = Bool.(seg_C_mask_minima)
    markers_C = label_components(seg_C_mask_minima)
    segment_C = watershed(distances_C, markers_C; mask=seg_C_mask_bool)
    labels_C = labels_map(segment_C)
    watershed_C = ifelse.(labels_C .== 0, 0, 1)

    ## Intersect the two watershed files
    watershed_intersect = watershed_B .* watershed_C

    return watershed_B, watershed_C, watershed_intersect
end