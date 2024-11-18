# Type of ref_image_band_2
__T__ = SubArray{
    Float64,
    2,
    Base.ReinterpretArray{Float64,3,RGB{Float64},Matrix{RGB{Float64}},true},
    Tuple{Int64,Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}}},
    false,
}

"""
    find_reflectance_peaks(reflectance_channel, possible_ice_threshold;)

Find histogram peaks in single channels of a reflectance image and return the second greatest peak. If needed, edges can be returned as the first object from `build_histogram`. Similarly, peak values can be returned as the second object from `findmaxima`.

# Arguments
- `reflectance_channel`: either band 2 or band 1 of false-color reflectance image
- `possible_ice_threshold`: threshold value used to identify ice if not found on first or second pass

"""
function find_reflectance_peaks(
    reflectance_channel::Union{__T__,Matrix{Float64}};
    possible_ice_threshold::Float64=Float64(75 / 255),
)::Int64
    reflectance_channel[reflectance_channel .< possible_ice_threshold] .= 0 #75 / 255
    _, counts = ImageContrastAdjustment.build_histogram(reflectance_channel)
    locs, _ = Peaks.findmaxima(counts)
    sort!(locs; rev=true)
    return locs[2] # second greatest peak
end

"""
    find_ice_labels(falsecolor_image, landmask; band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Locate the pixels of likely ice from false color reflectance image. Returns a binary mask with ice floes contrasted from background. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `falsecolor_image`: corrected reflectance false color image - bands [7,2,1]
- `landmask`: bitmatrix landmask for region of interest
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)

"""
function find_ice_labels(
    falsecolor_image::Matrix{RGB{Float64}},
    landmask::BitMatrix;
    band_7_threshold::Float64=Float64(5 / 255),
    band_2_threshold::Float64=Float64(230 / 255),
    band_1_threshold::Float64=Float64(240 / 255),
    band_7_threshold_relaxed::Float64=Float64(10 / 255),
    band_1_threshold_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)::Vector{Int64}

    ## Make ice masks
    cv = channelview(falsecolor_image)

    mask_ice_band_7 = @view(cv[1, :, :]) .< band_7_threshold #5 / 255
    mask_ice_band_2 = @view(cv[2, :, :]) .> band_2_threshold #230 / 255
    mask_ice_band_1 = @view(cv[3, :, :]) .> band_1_threshold #240 / 255
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice_labels = remove_landmask(landmask, ice)
    # @info "Done with masks" # to uncomment when logger is added

    ## Find likely ice floes
    if sum(abs.(ice_labels)) == 0
        mask_ice_band_7 = @view(cv[1, :, :]) .< band_7_threshold_relaxed #10 / 255
        mask_ice_band_1 = @view(cv[3, :, :]) .> band_1_threshold_relaxed #190 / 255
        ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
        ice_labels = remove_landmask(landmask, ice)
        if sum(abs.(ice_labels)) == 0
            ref_image_band_2 = @view(cv[2, :, :])
            ref_image_band_1 = @view(cv[3, :, :])
            band_2_peak = find_reflectance_peaks(ref_image_band_2, possible_ice_threshold = possible_ice_threshold)
            band_1_peak = find_reflectance_peaks(ref_image_band_1, possible_ice_threshold = possible_ice_threshold)
            mask_ice_band_2 = @view(cv[2, :, :]) .> band_2_peak / 255
            mask_ice_band_1 = @view(cv[3, :, :]) .> band_1_peak / 255
            ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
            ice_labels = remove_landmask(landmask, ice)
        end
    end
    # @info "Done with ice labels" # to uncomment when logger is added
    return ice_labels
end

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
    return missing
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
    factor,
)
Get the ice masks from the falsecolor image and morphological residue given a particualr tiling configuration.
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
    falsecolor_image,
    morph_residue,
    landmask::BitMatrix,
    tiles,
    binarize::Bool=true;
    band_7_threshold::T=5,
    band_2_threshold::T=230,
    band_1_threshold::T=240,
    band_7_threshold_relaxed::T=10,
    band_1_threshold_relaxed::T=190,
    possible_ice_threshold::T=75,
    factor::T=255,
) where {T<:Integer}

    # Make canvases
    ice_mask = BitMatrix(zeros(Bool, size(falsecolor_image)))
    binarized_tiling = zeros(Int, size(falsecolor_image))
    
    fc_landmasked = apply_landmask(falsecolor_image, landmask)

    for tile in tiles
        #  Conditionally update binarized_tiling as it's not used in some workflows
        if binarize
            binarized_tiling[tile...] .= imbinarize(morph_residue[tile...])
        end

        morph_residue_seglabels = kmeans_segmentation(Gray.(morph_residue[tile...] / 255))

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
    return (icemask=ice_mask, bin=binarized_tiling)
end
