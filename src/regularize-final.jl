using IceFloeTracker: unsharp_mask, to_uint8, hbreak, morph_fill

# TODO: Add tests for regularize_fill_holes, regularize_sharpening, get_final

"""
    regularize_fill_holes(img, local_maxima_mask, factor, segment_mask, L0mask)

Regularize `img` by:
    1. increasing the maxima of `img` by a factor of `factor`
    2. filtering `img` at positions where either `segment_mask` or `L0mask` are true
    3. filling holes

# Arguments
- `img`: The morphological residue image.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
- `L0mask`: zero-labeled pixels from watershed.
"""
function regularize_fill_holes(img, local_maxima_mask, segment_mask, L0mask, factor)
    new2 = to_uint8(img .+ local_maxima_mask .* factor)
    new2[segment_mask .|| L0mask] .= 0
    return IceFloeTracker.fill_holes(new2)
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
function regularize_sharpening(
    img, L0mask, local_maxima_mask, segment_mask, se, radius, amount, factor
)
    new3 = unsharp_mask(img, radius, amount, 255)
    new3[L0mask] .= 0
    new3 = IceFloeTracker.reconstruct(new3, se, "dilation", false)
    new3[segment_mask] .= 0
    return to_uint8(new3 + local_maxima_mask .* factor)
end

function _regularize(
    morph_residue, local_maxima_mask, segment_mask, L0mask, se; factor, radius, amount
)
    reg_fill_holes = regularize_fill_holes(
        morph_residue, local_maxima_mask, segment_mask, L0mask, factor[1]
    )
    reg_sharpened = regularize_sharpening(
        reg_fill_holes,
        L0mask,
        local_maxima_mask,
        segment_mask,
        se,
        radius,
        amount,
        factor[end],
    )
    return reg_sharpened
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
function get_final(
    # this function is used in preprocessing_tiling
    img,
    segment_mask,
    se_erosion,
    se_dilation,
    apply_segment_mask::Bool=true,
)
    _img = hbreak(img)

    # slow for big images
    _img .= morph_fill(_img)

    # TODO: decide on criteria for applying segment mask
    apply_segment_mask && (_img[segment_mask] .= false)

    # tends to fill more than matlabs imfill
    _img .= IceFloeTracker.fill_holes(_img)

    # marker image
    _img .= branch(_img)

    #= opening to remove noise while preserving shape/size
    Note the different structuring elements for erosion and dilation =#
    mask = sk_morphology.erosion(_img, se_erosion)
    mask .= sk_morphology.dilation(mask, se_dilation)

    # Restore shape of floes based on the cleaned up `mask`
    final = IceFloeTracker.mreconstruct(IceFloeTracker.dilate, _img, mask)
    return BitMatrix(final)
end
