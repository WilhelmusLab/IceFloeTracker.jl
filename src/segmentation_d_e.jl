"""
    segmentation_D()

Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl and segmentation_c.jl to further isolate ice floes, returning two binary segmentation masks.

# Arguments
- `segmentation_b_not_ice_mask`:
- `segmented_c`:


"""
function segmentation_D_E(
    segmentation_b_not_ice_mask::BitMatrix, segmented_c::BitMatrix;
)::Tuple{BitMatrix,BitMatrix}

    ## Watershed on segmentation_B intermediate
    feats_B = feature_transform(.!segmentation_b_not_ice_mask)
    seg_B = 1 .- (distance_transform(feats_B))
    seg_B_mask = ImageSegmentation.hmin_transform(seg_B, 2)
    seg_B_mask_minima = local_minima(seg_B_mask)
    seg_B_mask_minima[seg_B_mask_minima .> 0] .= 1
    seg_B_mask_bool = Bool.(seg_B_mask_minima)
    markers_B = label_components(seg_B_mask_minima)
    segment_B = watershed(seg_B, markers_B; mask=seg_B_mask_bool)
    labels_B = labels_map(segment_B)
    watershed_B = ifelse.(labels_B .== 0, 0, 1)

    ## Watershed on segmentation_C
    feats_C = feature_transform(.!segmented_c)
    seg_C = 1 .- (distance_transform(feats_C))
    seg_C_mask = ImageSegmentation.hmin_transform(seg_C, 2)
    seg_C_mask_minima = local_minima(seg_C_mask)
    seg_C_mask_minima[seg_C_mask_minima .> 0] .= 1
    seg_C_mask_bool = Bool.(seg_C_mask_minima)
    markers_C = label_components(seg_C_mask_minima)
    segment_C = watershed(seg_C, markers_C; mask=seg_C_mask_bool)
    labels_C = labels_map(segment_C)
    watershed_C = ifelse.(labels_C .== 0, 0, 1)

    return watershed_B, watershed_C
end
