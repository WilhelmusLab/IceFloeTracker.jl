using Images: build_histogram
using Peaks: findmaxima, peakproms!, peakwidths!
using DataFrames

# dmw: wrapper for the IceDetectionThreshold method. We should be able to 
# use this directly rather than piping the types.
# jgh: TODO: replace this function with the same call to binarize

# function get_ice_labels_mask(ref_img::Matrix{RGB{N0f8}}, thresholds)
#     mask =
#         binarize(
#             ref_img,
#             IceDetectionThresholdMODIS721(;
#                 band_7_max=thresholds[1],
#                 band_2_min=thresholds[2],
#                 band_1_min=thresholds[3],
#             ),
#         ) .|>
#         gray .|>
#         Bool
#     @debug "Found $(sum(mask)) ice pixels"
#     return mask
# end

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
    
    f = IceDetectionFirstNonZeroAlgorithm([
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_threshold,
            band_2_min=band_2_threshold,
            band_1_min=band_1_threshold
        ),
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_threshold_relaxed,
            band_2_min=band_2_threshold,
            band_1_min=band_1_threshold_relaxed,
        ),
        IceDetectionBrightnessPeaksMODIS721(;
            band_7_max=band_7_threshold,
            possible_ice_threshold=possible_ice_threshold
        ),
        IceDetectionThresholdMODIS721(;
            band_7_max=1.,
            band_2_min=band_2_threshold,
            band_1_min=0.,
        )
        ])

    ice_labels = binarize(falsecolor_img, f)
    isempty(ice_labels) && return -1
    sum(ice_labels) == 0 && return -1
    return StatsBase.mode(image_indexmap[mask])
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
- Binary image with likely sea ice floes = 1.
"""
function get_ice_masks( #tbd: rename to kmeans_binarization?
    falsecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    morph_residue::AbstractArray{<:AbstractGray},
    landmask::AbstractArray{<:Bool},
    tiles::AbstractMatrix{Tuple{UnitRange{Int64},UnitRange{Int64}}},
    binarize::Bool=true, # remove if the get_nlabel works
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
