using IceFloeTracker: unsharp_mask, to_uint8, reconstruct, hbreak, to_uint8, morph_fill
using IceFloeTracker.MorphSE: dilate, erode, fill_holes

# TODO: Add tests for regularize_fill_holes, regularize_sharpening, get_final

"""
    regularize_fill_holes(img, local_maxima_mask, factor, segment_mask, L0mask)

Regularize `img` by:
    1. increasing the maxima of `img` by a factor of `factor`
    2. filtering `img` at positions where `segment_mask` and `L0mask` are true
    3. filling holes

# Arguments
- `img`: The morphological residue image.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
- `L0mask`: zero-labeled pixels from watershed.
"""
function regularize_fill_holes(img, local_maxima_mask, factor, segment_mask, L0mask)
    new2 = to_uint8(img .+ local_maxima_mask .* factor)
    new2[segment_mask .|| L0mask] .= 0
    return fill_holes(new2)
end

"""
    regularize_sharpening(img, L0mask, radius, amount, local_maxima_mask, factor, segment_mask, se)

Regularize `img` via sharpening, filtering, reconstruction, and maxima elevating.

# Arguments
- `img`: The input image.
- `L0mask`: zero-labeled pixels from watershed.
- `radius`: The radius of the unsharp mask.
- `amount`: The amount of unsharp mask.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
"""
function regularize_sharpening(img, L0mask, radius, amount, local_maxima_mask, factor, segment_mask, se)
    new3 = unsharp_mask(img, radius, amount, 255)
    new3[L0mask] .= 0
    new3 = reconstruct(new3, se, "dilation", false)
    new3[segment_mask] .= 0
    return to_uint8(new3 + local_maxima_mask .* factor)
end

"""
    get_final(img, label, segment_mask, se_erosion, se_dilation)

Final processing following the tiling workflow.

# Arguments
- `img`: The input image.
- `label`: Mode of most common label in the find_ice_labels workflow.
- `segment_mask`: The segment mask.
- `se_erosion`: structuring element for erosion.
- `se_dilation`: structuring element for dilation.
- `apply_segment_mask=true`: Whether to filter `img` the segment mask.

"""
function get_final(img, label, segment_mask, se_erosion, se_dilation)
    img = hbreak(img)

    # slow for big images
    img = morph_fill(img)

    # only works for label 1, whose value tends to be arbirary.
    # Added for consistency with MASTER.m. CP
    if label == 1
        img[segment_mask] .= false
    end

    # tends to fill more than matlabs imfill
    img = fill_holes(img)

    marker = branch(img)

    mask = erode(marker, se_erosion)
    mask = dilate(mask, se_dilation)

    # Added for consistency with MASTER.m. CP
    if label == 1
        mask[1] = false
    end

    final = sk_morphology.reconstruction(marker, mask)
    return final
end
