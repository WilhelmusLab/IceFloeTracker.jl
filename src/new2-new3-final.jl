using IceFloeTracker: MorphSE, unsharp_mask, to_uint8, reconstruct

# TODO: choose better names

# TODO: Add tests for get_new2, get_new3, get_final

"""
    getnew2(morph_residue, local_maxima_mask, factor, segment_mask, L0mask)
Calculate the new image `new2` from the input image `morph_residue`.
# Arguments
- `morph_residue`: The morphological residue image.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
- `L0mask`: zero-labeled pixels from watershed.
"""
function get_new2(morph_residue, local_maxima_mask, factor, segment_mask, L0mask)
    new2 = to_uint8(morph_residue .+ local_maxima_mask .* factor)
    new2[segment_mask .|| L0mask] .= 0
    return MorphSE.fill_holes(new2)
end

"""
    get_new3(new2, L0mask, radius, amount, local_maxima_mask, factor, segment_mask)
Calculate the new image `new3` from the input image `new2`.
# Arguments
- `img`: The input image.
- `L0mask`: zero-labeled pixels from watershed.
- `radius`: The radius of the unsharp mask.
- `amount`: The amount of unsharp mask.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
"""
function get_new3(img, L0mask, radius, amount, local_maxima_mask, factor, segment_mask)
    new3 = unsharp_mask(img, radius, amount, 255)
    new3[L0mask] .= 0
    new3 = reconstruct(new3, se, "dilation", false)
    new3[segment_mask] .= 0
    return to_uint8(new3 + local_maxima_mask .* factor)
end

"""
    get_final(img, label, segment_mask, se_erosion, se_dilation)

Final processing following the tiling workflow.
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
    img = IceFloeTracker.MorphSE.fill_holes(img)

    marker = branch(img)

    mask = MorphSE.erode(marker, se_erosion)
    mask = MorphSE.dilate(mask, se_dilation)

    # Added for consistency with MASTER.m. CP
    if label == 1
        mask[1] = false
    end

    final = sk_morphology.reconstruction(marker, mask)
    return final
end
