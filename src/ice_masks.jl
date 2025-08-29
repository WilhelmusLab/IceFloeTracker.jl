using Images: build_histogram
using Peaks: findmaxima, peakproms!, peakwidths!
using DataFrames

# """
# Given the edges and counts from build_histogram, identify local maxima and return the location of the
# largest local maximum that is bright enough that it is possibly sea ice. Locations are determined by 
# the edges, which by default are the left bin edges. Note also that peaks defaults to the left side of
# plateaus.
# """
# function get_ice_peaks(edges, counts; possible_ice_threshold::Float64=0.30, minimum_prominence::Float64=0.05, window::Int64=3)
#     counts = counts[1:end]
#     counts = counts ./ sum(counts[edges .> possible_ice_threshold])
#     pks = findmaxima(counts, window) |> peakproms! |> peakwidths!
#     pks_df = DataFrame(pks[Not(:data)])
#     pks_df = sort(pks_df, :proms, rev=true)
#     maximum(pks_df.proms) < minimum_prominence && return Inf
#     return edges[pks_df[argmax(pks_df.proms), :indices]] 
# end

# dmw: wrapper for the IceDetectionThreshold method. We should be able to 
# use this directly rather than piping the types.
# jgh: TODO: replace this function with the same call to binarize
function get_ice_labels_mask(ref_img::Matrix{RGB{N0f8}}, thresholds)
    mask =
        binarize(
            ref_img,
            IceDetectionThresholdMODIS721(;
                band_7_max=thresholds[1],
                band_2_min=thresholds[2],
                band_1_min=thresholds[3],
            ),
        ) .|>
        gray .|>
        Bool
    @debug "Found $(sum(mask)) ice pixels"
    return mask
end

function get_nlabel(
    falsecolor_img,
    segmented_image_indexmap;
    band_7_threshold::T=5/255,
    band_2_threshold::T=230/255,
    band_1_threshold::T=240/255,
    band_7_threshold_relaxed::T=10/255,
    band_1_threshold_relaxed::T=190/255,
    possible_ice_threshold::T=75/255,
) where {T<:Float64}
    _getnlabel(image_indexmap, mask) = begin
        isempty(mask) && return -1
        sum(mask) == 0 && return -1
        StatsBase.mode(image_indexmap[mask])
    end

    # Initial threshold set identifies bright ice
    thresholds = (band_7_threshold, band_2_threshold, band_1_threshold)
    ice_labels_mask = get_ice_labels_mask(falsecolor_img, thresholds)
    sum(ice_labels_mask) > 0 && return _getnlabel(segmented_image_indexmap, ice_labels_mask)
    @debug "Trying first relaxation."

    # First relaxation allows slightly grayer ice, potentially with some cloud
    thresholds = (band_7_threshold_relaxed, band_2_threshold, band_1_threshold_relaxed)
    ice_labels_mask = get_ice_labels_mask(falsecolor_img, thresholds)
    sum(ice_labels_mask) > 0 && return _getnlabel(segmented_image_indexmap, ice_labels_mask)

    @debug "Trying second/third relaxation."
    # The second and third relaxation are handled in a separate function.
    return get_nlabel_relaxation(
        falsecolor_img,
        segmented_image_indexmap,
        possible_ice_threshold, # set this with data
        band_7_threshold_relaxed, # set this with data
        band_2_threshold,
    )
end

function get_nlabel_relaxation(
    falsecolor_img,
    segmented_image_indexmap,
    possible_ice_threshold,
    band_7_threshold_relaxed,
    band_2_threshold,
)
    _getnlabel(image_indexmap, mask) = begin
        isempty(mask) && return -1
        sum(mask) == 0 && return -1
        StatsBase.mode(image_indexmap[mask])
    end

    # filter b/c channels (landmasked channels 2 and 3) and compute peaks
    band_2 = green.(falsecolor_img)
    band_1 = blue.(falsecolor_img)
    # band_2, band_1 = [float64.(channelview(falsecolor_img)[i, :, :]) for i in 2:3]

    # temporary fix in case of integer thresholds
    p = possible_ice_threshold > 1 ? possible_ice_threshold / 255 : possible_ice_threshold
    band_2[band_2 .< p] .= 0
    band_1[band_1 .< p] .= 0

    pks_band_2 = get_ice_peaks(build_histogram(band_2, 64; minval=0, maxval=1)...)
    pks_band_1 = get_ice_peaks(build_histogram(band_1, 64; minval=0, maxval=1)...)

    # pks_band_2, pks_band_1 = get_image_peaks.([band_2, band_1])

    # return early if no peaks are found
    isinf(pks_band_2) || isinf(pks_band_1) && return -1
    # !all(length.([pks_band_2.locs, pks_band_1.locs]) .> 2) && return -1

    # relaxed_thresholds = [band_7_threshold_relaxed, pks_band_2.locs[2], pks_band_1.locs[2]]
    relaxed_thresholds = [band_7_threshold_relaxed, pks_band_2, pks_band_1]
    ice_labels_mask = get_ice_labels_mask(falsecolor_img, relaxed_thresholds)

    sum(ice_labels_mask) > 0 && return _getnlabel(segmented_image_indexmap, ice_labels_mask)

    # Final relaxation
    ice_labels_mask = band_2 .> band_2_threshold
    return _getnlabel(segmented_image_indexmap, ice_labels_mask)
end

# dmw: split into the k-means and binarization methods, since they operate on different principles.
# remove the "factor" argument, since it can be inferred from the image type.
# move this into a segmentation algorithms file
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
        k
    )

Identifies potential sea ice floes using two methods: selection of a relevant k-means cluster and application of
adaptive threshold binarization. For the k-means section, a series of thresholds on band 7, 2, and 1 reflectance 
are applied in order to find the cluster containing bright sea ice pixels.

# Arguments
- `falsecolor_image`: MODIS False Color Bands 7-2-1.
- `morph_residue`: Grayscale sharpened and equalized image from preprocessing workflow.
- `landmask`: Binary landmask. 
- `tiles`: Iterable with tile divisions.
- `binarize::Bool=true`: Whether to binarize the tiling.
- `band_7_threshold=5/255`: The threshold for band 7.
- `band_2_threshold=230/255`: The threshold for band 2.
- `band_1_threshold=240/255`: The threshold for band 1.
- `band_7_threshold_relaxed=10`: The relaxed threshold for band 7.
- `band_1_threshold_relaxed=190`: The relaxed threshold for band 1.
- `possible_ice_threshold=75/255`: The threshold for possible ice.
- `k=4`: The number of clusters to use for k-means segmentation.

# Returns
- A named tuple `(icemask, bin)` where:
  - `icemask`: The ice mask.
  - `bin`: The binarized tiling.
"""
function get_ice_masks(
    falsecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    morph_residue::AbstractArray{<:AbstractGray},
    landmask::AbstractArray{<:Bool},
    tiles::AbstractMatrix{Tuple{UnitRange{Int64},UnitRange{Int64}}},
    # TODO: don't shadow "binarize" function (dmw: remove binarization from this function)
    binarize::Bool=true;
    band_7_threshold::Float64=5/255,
    band_2_threshold::Float64=230/255,
    band_1_threshold::Float64=240/255,
    band_7_threshold_relaxed::Float64=10/255,
    band_1_threshold_relaxed::Float64=190/255,
    possible_ice_threshold::Float64=75/255,
    k=4
) 

    # Make canvases
    sz = size(falsecolor_image)
    ice_mask = BitMatrix(zeros(Bool, sz))
    binarized_tiling = zeros(Int, sz)

    fc_landmasked = apply_landmask(falsecolor_image, .!landmask) # will need to flip once the landmask is the right style

    # Threads.@threads
    for tile in tiles
        @debug "Processing tile: $tile"
        mrt = morph_residue[tile...]
        segmented_image_indexmap = kmeans_segmentation(mrt; k=k)

        # TODO: handle case where get_nlabel returns missing
        floes_label = get_nlabel(
            fc_landmasked[tile...],
            segmented_image_indexmap,           
            band_7_threshold=band_7_threshold,
            band_2_threshold=band_2_threshold,
            band_1_threshold=band_1_threshold,
            band_7_threshold_relaxed=band_7_threshold_relaxed,
            band_1_threshold_relaxed=band_1_threshold_relaxed,
            possible_ice_threshold=possible_ice_threshold,
        )

        ice_mask[tile...] .= (segmented_image_indexmap .== floes_label)

        #  Conditionally update binarized_tiling as its not used in some workflows
        binarize && (binarized_tiling[tile...] .= imbinarize(mrt))
    end

    return (icemask=ice_mask, bin=binarized_tiling .> 0)
end
