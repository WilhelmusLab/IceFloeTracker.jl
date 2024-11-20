function get_image_peaks(arr, imgtype="uint8")
    _, heights = imhist(arr, imgtype)

    locs, heights, _ = Peaks.findmaxima(heights)

    # TODO: make this conditional on input args
    order = sortperm(heights; rev=true)
    locs, heights = locs[order], heights[order]

    return (locs=locs, heights=heights)
end

function get_ice_labels_mask(ref_img::Matrix{RGB{N0f8}}, thresholds, factor=1)
    cv = channelview(ref_img)
    cv = [float64.(cv[i, :, :]) .* factor for i in 1:3]
    mask_ice_band_7 = cv[1] .< thresholds[1]
    mask_ice_band_2 = cv[2] .> thresholds[2]
    mask_ice_band_1 = cv[3] .> thresholds[3]
    mask = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    @debug "Found $(sum(mask)) ice pixels"
    return mask
end

function get_nlabel(
    ref_img,
    morph_residue_labels,
    factor;
    band_7_threshold::T=5,
    band_2_threshold::T=230,
    band_1_threshold::T=240,
    band_7_threshold_relaxed::T=10,
    band_1_threshold_relaxed::T=190,
    possible_ice_threshold::T=75,
) where {T<:Integer}
    _getnlabel(morphr, mask) = StatsBase.mode(morphr[mask])

    # Initial attempt to get ice labels
    thresholds = (band_7_threshold, band_2_threshold, band_1_threshold)
    ice_labels_mask = get_ice_labels_mask(ref_img, thresholds, 255)
    sum(ice_labels_mask) > 1 && return _getnlabel(morph_residue_labels, ice_labels_mask)

    # First relaxation
    thresholds = (band_7_threshold_relaxed, band_2_threshold, band_1_threshold_relaxed)
    ice_labels_mask = get_ice_labels_mask(ref_img, thresholds, 255)
    sum(ice_labels_mask) > 0 && return _getnlabel(morph_residue_labels, ice_labels_mask)

    # Second/Third relaxation
    return get_nlabel_relaxation(
        ref_img,
        morph_residue_labels,
        factor,
        possible_ice_threshold,
        band_7_threshold_relaxed,
        band_2_threshold,
    )
end

function get_nlabel_relaxation(
    ref_img,
    morph_residue_labels,
    factor,
    possible_ice_threshold,
    band_7_threshold_relaxed,
    band_2_threshold,
)
    # filter b/c channels (landmasked channels 2 and 3) and compute peaks
    b, c = [float64.(channelview(ref_img)[i, :, :]) .* factor for i in 2:3]
    b[b .< possible_ice_threshold] .= 0
    c[c .< possible_ice_threshold] .= 0
    pksb, pksc = get_image_peaks.([b, c])

    # return early if no peaks are found
    !all(length.([pksb.locs, pksc.locs]) .> 2) && return 1

    relaxed_thresholds = [band_7_threshold_relaxed, pksb.locs[2], pksc.locs[2]]
    ice_labels = get_ice_labels_mask(ref_img, relaxed_thresholds, factor)

    sum(ice_labels) > 0 && return StatsBase.mode(morph_residue_labels[ice_labels])

    # Final relaxation
    mask_b = b .> band_2_threshold
    sum(mask_b) > 0 && return StatsBase.mode(morph_residue_labels[mask_b])

    # No mode found
    return 1
end

"""
    get_ice_masks(
        falsecolor_image,
        morph_residue,
        landmask,
        tiles,
        binarize;
        band_7_threshold,
        band_2_threshold,
        band_1_threshold,
        band_7_threshold_relaxed,
        band_1_threshold_relaxed,
        possible_ice_threshold,
        factor
    )

Get the ice masks from the falsecolor image and morphological residue given a particular tiling configuration.

# Arguments
- `falsecolor_image`: The falsecolor image.
- `morph_residue`: The morphological residue image.
- `landmask`: The landmask.
- `tiles`: The tiles.
- `binarize::Bool=true`: Whether to binarize the tiling.
- `band_7_threshold=5`: The threshold for band 7.
- `band_2_threshold=230`: The threshold for band 2.
- `band_1_threshold=240`: The threshold for band 1.
- `band_7_threshold_relaxed=10`: The relaxed threshold for band 7.
- `band_1_threshold_relaxed=190`: The relaxed threshold for band 1.
- `possible_ice_threshold=75`: The threshold for possible ice.
- `factor=255`: normalization factor to convert images to uint8.
# Returns
- A named tuple `(icemask, bin)` where:
  - `icemask`: The ice mask.
  - `bin`: The binarized tiling.
"""
function get_ice_masks(
    falsecolor_image::Matrix{RGB{N0f8}},
    morph_residue::Matrix{<:Integer},
    landmask::BitMatrix,
    tiles::S,
    binarize::Bool=true;
    band_7_threshold::T=5,
    band_2_threshold::T=230,
    band_1_threshold::T=240,
    band_7_threshold_relaxed::T=10,
    band_1_threshold_relaxed::T=190,
    possible_ice_threshold::T=75,
    factor::T=255,
) where {T<:Integer,S<:AbstractMatrix{Tuple{UnitRange{Int64},UnitRange{Int64}}}}

    # Make canvases
    sz = size(falsecolor_image)
    ice_mask = BitMatrix(zeros(Bool, sz))
    binarized_tiling = zeros(Int, sz)

    fc_landmasked = apply_landmask(falsecolor_image, landmask)

    # Threads.@threads
    for tile in tiles
        #  Conditionally update binarized_tiling as its not used in some workflows
        if binarize
            binarized_tiling[tile...] .= imbinarize(morph_residue[tile...])
        end

        morph_residue_seglabels = kmeans_segmentation(Gray.(morph_residue[tile...] / 255))

        # TODO: handle case where get_nlabel returns missing
        floes_label = get_nlabel(
            fc_landmasked[tile...],
            morph_residue_seglabels,
            factor;
            band_7_threshold=band_7_threshold,
            band_2_threshold=band_2_threshold,
            band_1_threshold=band_1_threshold,
            band_7_threshold_relaxed=band_7_threshold_relaxed,
            band_1_threshold_relaxed=band_1_threshold_relaxed,
            possible_ice_threshold=possible_ice_threshold,
        )

        ice_mask[tile...] .= (morph_residue_seglabels .== floes_label)
    end

    return (icemask=ice_mask, bin=binarized_tiling .> 0)
end
