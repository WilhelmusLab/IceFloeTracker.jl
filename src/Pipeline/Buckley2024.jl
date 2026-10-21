module Buckley2024

using Images
using DataFrames

import ..Pipeline: 
    IceFloeSegmentationAlgorithm

import ..Segmentation:
    erosion_expansion,
    view_seg,
    view_seg_random

"""
    Buckley2024.Segment()

Implementation of the processing method for Buckley et al. 2024.

# Parameters
- `coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))`: Structuring element for the `create_landmask` function
- `cloud_mask_algorithm = Watkins2025CloudMask()`: Cloud mask algorithm
- `preprocessing_algorithm = Preprocess()`: Function to sharpen and equalize the truecolor image
- `tile_size_pixels=1200`: Nominal tile size in pixels
- `min_tile_ice_pixel_count=300`: Smallest number of required sea ice pixels in tile
- `preliminary_ice_mask = IceDetectionBrightnessMidpoint(minimum_reflectance=0.3)`: Function to use to identify likely ice pixels for filtering.
- `kmeans_params = (k=4, maxiter=50, random_seed=45)`: Parameters for `kmeans_binarization`
- `cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
    band_7_max=0.1,
    possible_ice_threshold=0.3,
    join_method="union",
    minimum_prominence=0.01)`: Function to use to select a k-means cluster in the `kmeans_binarization` workflow
- `clean_binary_floes_params`: Parameters for the preliminary binary image cleanup
- `floe_splitting_params`: Parameters for the `dist_morph_split` floe splitting algorithm
- `floe_filtering_params`: Parameters for post-segmentation cleanup
"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} =
        coastal_buffer_structuring_element
    cloud_mask_algorithm = cloud_mask_algorithm
    preprocessing_algorithm = preprocessing_algorithm
    preliminary_ice_mask = preliminary_ice_mask
    kmeans_params = kmeans_params
    adaptive_params = adaptive_params
    floe_splitting_params = floe_splitting_params
    floe_filtering_params = floe_filtering_params
end

function (s::Segment)(
    truecolor::T₁,
    falsecolor::T₂,
    landmask::T₃,
    coastal_buffer_mask::T₄;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {
    T₁<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₂<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₃<:AbstractMatrix{<:Union{Bool,Gray{Bool}}},
    T₄<:AbstractMatrix{<:Union{Bool,Gray{Bool}}},
}
    # Move these conversions down through the function as each step gets support for
    # the full range of image formats
    truecolor_image = float64.(RGB.(truecolor))
    falsecolor_image = float64.(RGB.(falsecolor)) # Only used for cloud mask
    landmask = landmask .> 0 # make sure it's a bitmatrix, in case it's passed as Gray
    apply_landmask!(truecolor_image, landmask)
    apply_landmask!(falsecolor_image, landmask)

    cloud_mask = create_cloudmask(falsecolor_image, s.cloud_mask_algorithm)
    tc_masked = apply_landmask(truecolor, cloud_mask)

    # Ice mask from red channel of truecolor image
    gray_img = Gray.(red.(tc_masked))
    ice_mask = binarize(gray_img, AdaptiveThreshold(; window_size=399, percentage=0))
    # Also get the min/max thresholds for the image 
    apply_landmask!(ice_mask, landmask .|| cloud_mask)

    @info "Splitting floes"
    final_floes = erosion_expansion(
        gray_img, 
        ice_mask; 
        itmin=1, # minimum number of iterations
        itmax=8, # maximum number of iterations
        step=1, # spacing between iterations
        strel=strel_diamond((3,3)),
        minimum_area=16,
    )

    # filter steps: remove dark floes, remove small floes
    # TODO: remove segments if segment mean is too low.
    remove_small_segments!(final_floes, s.floe_filtering_params.min_floe_size)
    remove_large_segments!(final_floes, s.floe_filtering_params.max_floe_size)

    # Re-label so there are no missing numbers in the component list
    final_floes .= label_components(final_floes)

    # Return the original truecolor image, segmented
    segments_tc = SegmentedImage(truecolor_image, final_floes)
    segments_fc = SegmentedImage(falsecolor_image, final_floes)

    if !isnothing(intermediate_results_callback)
        colorview_random = view_seg_random(segments_tc)
        segment_mean_truecolor=n0f8.(segment_mean_map(segments_tc))
        segment_mean_falsecolor=n0f8.(segment_mean_map(segments_fc))
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            coastal_buffer_mask=Gray.(coastal_buffer_mask),
            cloud_mask=Gray.(cloud_mask),
            ice_mask=Gray.(prelim_ice_mask),
            final_floes=colorview_random,
            labels_map=final_floes,
            segment_mean_falsecolor=segment_mean_falsecolor,
            segment_mean_truecolor=segment_mean_truecolor,
        )
    end
    return segments_tc
end



end