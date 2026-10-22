module Buckley2024

using Images
using DataFrames

import ..Pipeline: 
    IceFloeSegmentationAlgorithm
import ..Preprocessing:
    create_landmask,
    create_cloudmask,
    apply_landmask,
    apply_landmask!,
    Watkins2026CloudMask

import ..Segmentation:
    erosion_expansion,
    remove_large_segments!,
    remove_small_segments!,
    view_seg,
    view_seg_random

# Parameters
adaptive_params = (; window_size=399, percentage=0)
coastal_buffer_structuring_element = strel_box((25, 25))
cloud_mask_algorithm = Watkins2026CloudMask()
floe_splitting_params = (
    min_floe_area=16,
    iter_max=8,
    iter_min=1,
    iter_step=1,
    strel=strel_diamond((3,3)),
    )
floe_filtering_params = (
    min_floe_size=16,
    max_floe_size=90e3, # TODO: Check if Ellen used a maximum floe size
    min_band_1_reflectance=0.58,
    min_circularity=0., # TODO: Check if Ellen used this threshold
)
"""
    Buckley2024.Segment()

Implementation of the processing method for Buckley et al. 2024.

# Parameters
- `coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))`: Structuring element for the `create_landmask` function
- `cloud_mask_algorithm = Watkins2026CloudMask()`: Cloud mask algorithm
- `adaptive_params`: parameters to send to the AdaptiveThreshold binarization method
- `floe_splitting_params`: Parameters for the `erosion_expansion` floe splitting algorithm
- `floe_filtering_params`: Parameters for post-segmentation cleanup
"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} =
        coastal_buffer_structuring_element
    cloud_mask_algorithm = cloud_mask_algorithm
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
    ice_mask = binarize(gray_img, AdaptiveThreshold(s.adaptive_params...))
    # Also get the min/max thresholds for the image 
    apply_landmask!(ice_mask, landmask .|| cloud_mask)

    @info "Splitting floes"
    final_floes = erosion_expansion(
        gray_img, 
        ice_mask; 
        itmin=s.floe_splitting_params.iter_min,
        itmax=s.floe_splitting_params.iter_max,
        step=s.floe_splitting_params.iter_step, 
        strel=s.floe_splitting_params.strel,
        minimum_area=s.floe_splitting_params.min_floe_area,
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