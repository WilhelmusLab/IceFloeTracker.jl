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
  - `label`: Most frequent label in the ice mask.
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

    Threads.@threads for tile in tiles
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