"""

This module documents the simplified IFT pipeline calibrated and validated in the paper
"Automatic segmentation and tracking of sea ice floes using Ice Floe Tracker algorithm" by
Watkins et al., submitted to IEEE Transactions in Geoscience and Remote Sensing (TGRS).

The pipeline contains distinct struct/functor pairs and uses both semantic and object-based 
analysis techniques. It introduces Preprocess() and Classify() functors, extending the Segment()
and Track() functors in LopezAcosta2019 and LopezAcosta2019Tiling.

"""

module TGRS2026

import ..Filtering:
    ContrastLimitedAdaptiveHistogramEqualization

import ..ImageUtils: 
    get_tiles, 
    apply_landmask,
    apply_landmask! # TODO: Test using the "masker" approach

import ..Preprocessing:
    Watkins2026CloudMask

import ..Segmentation:
    dist_morph_split,
    get_relevant_set,
    kmeans_segmentation,
    IceDetectionBrightnessMidpoint,
    regionprops_table,
    remove_small_segments!,
    remove_large_segments!,
    stitch_clusters,
    view_seg,
    view_seg_random

import ..Tracking:
    ChainedFilterFunction,
    DistanceThresholdFilter,
    euclidean_distance,
    FloeTracker,
    LogLogQuadraticTimeDistanceFunction,
    MinimumWeightMatchingFunction,
    PiecewiseLinearThresholdFunction,
    RelativeErrorThresholdFilter,
    ShapeDifferenceThresholdFilter,
    PsiSCorrelationThresholdFilter

import ..Pipeline: IceFloeSegmentationAlgorithm

abstract type IceFloePreprocessingAlgorithm end
abstract type IceFloeClassificationAlgorithm end

"""
   Preprocess(
        adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=1)
    )
    Preprocess()(img, mask)

    Converts input image to grayscale, then preprocesses by applying contrast limited adaptive histogram
    equalization. The mask may include the land mask, coastal buffer, or a domain

"""
@kwdef struct Preprocess <: IceFloePreprocessingAlgorithm
    histogram_algorithm = ContrastLimitedAdaptiveHistogramEqualization
    histogram_params = (nbins=256, rblocks=4, cblocks=4, clip=1)
end

function (p::Preprocess)(
    image::AbstractArray{<:Union{AbstractGray, TransparentGray, AbstractRGB,TransparentRGB}}, landmask
)
    # Cast to grayscale first to save compute time
    proc_img = Gray.(image)
    apply_landmask!(proc_img, landmask)

    adjust_histogram!(
        proc_img,
        p.histogram_algorithm(;
            p.histogram_params...
        ),
    )

    # Re-apply mask so histogram adjustment doesn't bleed into land
    apply_landmask!(proc_img, landmask)
    return proc_img
end

"""
   Classify(
        τ₁=0.1,
        τ₂=0.2,
        τ₇=0.2,
        key=Dict("land"=>0, "water"=>1, "ice"=>2, "cloud"=>3)
    )
    Classify()(false_color_image, mask)

Classifies an image into land, water, ice, and cloud using the Watkins2026 cloud mask
and the IceDetectionBrightnessMidpoint algorithm. The parameter τ₁ is the brightness 
minimum for the ice detection algorithm, while the τ₂ and τ₇ parameters are used in the 
cloud mask algorithm. The `key` specifies the integers used to encode the classification
for the returned label map.

"""
@kwdef struct Classify <: IceFloeClassificationAlgorithm
        τ₁=0.1
        τ₂=0.2
        τ₇=0.2
        key=Dict("land"=>0, "water"=>1, "ice"=>2, "cloud"=>3)
end

function (c::IceFloeClassificationAlgorithm)(false_color_image, land_mask)::Matrix{Int64}
    cloud_mask_algorithm=Watkins2026CloudMask(band_2_threshold=c.τ₂, band_7_threshold=c.τ₇)
    ice_mask_algorithm=IceDetectionBrightnessMidpoint(; minimum_reflectance=c.τ₁)
    fc_masked = apply_landmask(false_color_image, land_mask)
    clouds = cloud_mask_algorithm(fc_masked)
    ice = Gray.(blue.(apply_landmask(fc_masked, clouds))) .|> ice_mask_algorithm
    
    classified_image = ones(Int64, size(false_color_image)) .* c.key["water"]
    classified_image[land_mask .> 0] .= c.key["land"]
    classified_image[ice .> 0] .= c.key["ice"]
    classified_image[clouds .> 0] .= c.key["cloud"]
    
    return classified_image
end

"""
    Segment()

Function to identify pixel types and detect ice floes from MODIS true color and false color imagery.

"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_disk(5)
    tile_settings = (; rblocks=1, cblocks=1)
    min_ice_pixel_count = 300
    preprocessing_algorithm = Preprocess()
    classification_algorithm = Classify()
    classification_key = Dict("land"=>0, "water"=>1, "ice"=>2, "cloud"=>3) # TODO: way to share this? (maybe just classification_algorithm.key?)
    kmeans_params = kmeans_params # DetectFloes() functor, which has the clean, split, and filter on board?
    cleanup_binary_params = clean_binary_floes
    floe_splitting_algorithm = dist_morph_split
    floe_filtering_params = filter_floes
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
    falsecolor_image = float64.(RGB.(falsecolor))
    landmask = landmask .> 0 # make sure it's a bitmatrix, in case it's passed as Gray

    # n, m = size(truecolor_image) TODO: warn if the tile size ends up smaller than 400 px
    tiles = get_tiles(truecolor; p.tile_settings...)
    
    @info "Preprocess"
    preproc_gray = p.preprocessing_algorithm(truecolor_image, landmask)

    @info "Classify"
    classified_image = p.classification_algorithm(falsecolor_image, landmask)

    # Then check for sufficient possible sea ice pixels
    ice_mask = classified_image .== p.classification_algorithm.key["ice"]
    filtered_tiles = filter(
        t -> sum(ice_mask[t...]) > s.min_ice_pixel_count, filtered_tiles
    );

    @info "Detect Floes"
    # We use the cloud mask in finding the bright floes - the bright floe cluster can't be cloud -
    # and allow the k-means cluster to overlap with the cloud mask by using the preproc gray with
    # only the landmask applied to it. Not applying the cloudmask to the kmeans result, though, means
    # we need to be careful about the clouds.
    kmeans_result = kmeans_binarization(
        preproc_gray, fc_masked, filtered_tiles; s.kmeans_params...
    )
    adaptive_result = binarize(preproc_gray, AdaptiveThreshold(; s.adaptive_params...)) .> 0

    # AdaptiveThreshold often has noise in large blank areas
    apply_landmask!(adaptive_result, landmask)

    # We also don't want to include artificially brightened regions, so
    # we mask things that have already been classified as water.
    apply_landmask!(adaptive_result, .!(prelim_ice_mask .|| cloud_mask))

    @info "Splitting floes"
    clean_split_label =
        r -> dist_morph_split(
            clean_binary_floes(r, prelim_ice_mask, cloud_mask; s.cleanup_binary_params...);
            s.floe_splitting_params...,
        )

    kmeans_split_floes = clean_split_label(kmeans_result)
    adaptive_split_floes = clean_split_label(adaptive_result)

    # TBD: Filter floes based on the edge properties, colors

    @info "Filtering floes"
    filter_floes!(
        kmeans_split_floes,
        coastal_buffer_mask,
        cloud_mask,
        falsecolor_image;
        s.floe_filtering_params...,
    )
    filter_floes!(
        adaptive_split_floes,
        coastal_buffer_mask,
        cloud_mask,
        falsecolor_image;
        s.floe_filtering_params...,
    )

    @info "Joining segmentation results"
    final_floes = merge_floes(kmeans_split_floes, adaptive_split_floes, preproc_gray)

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
            preprocessed=preproc_gray,
            binarized=kmeans_result .> 0,
            final_floes=colorview_random,
            labels_map=final_floes,
            segment_mean_falsecolor=segment_mean_falsecolor,
            segment_mean_truecolor=segment_mean_truecolor,
        )
    end
    return segments_tc
end




#### Tracker parameters ####

const max_travel_distance_filter = DistanceThresholdFilter(;
    threshold_function=LogLogQuadraticTimeDistanceFunction()
)

const area_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:area,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.43, maximum_value=0.17
    ),
)

const convex_area_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:convex_area,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.44, maximum_value=0.25
    ),
)

const major_axis_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:major_axis_length,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.27, maximum_value=0.13
    ),
)

const minor_axis_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:minor_axis_length,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.28, maximum_value=0.1
    ),
)

const shape_difference_filter = ShapeDifferenceThresholdFilter(;
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.47, maximum_value=0.31
    ),
)

const psi_s_correlation_filter = PsiSCorrelationThresholdFilter(;
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.86, maximum_value=0.96
    ),
)

const FSFilterFunctions = [
    max_travel_distance_filter,
    area_relative_error_filter,
    convex_area_relative_error_filter,
    major_axis_relative_error_filter,
    minor_axis_relative_error_filter,
    shape_difference_filter,
    psi_s_correlation_filter,
]

const FSMatchingColumns = [
            :scaled_distance,
            :relative_error_area,
            :relative_error_convex_area,
            :relative_error_major_axis_length,
            :relative_error_minor_axis_length,
            :psi_s_correlation_score,
            :scaled_shape_difference,
        ]
"""
    Track()

Track shapes across images using the LogLogQuadratic distance filter, the ChainedFilterFunction,
and the MinimumWeightMatchingFunction.

"""
function Track(
    filter_function=ChainedFilterFunction(; filters=FSFilterFunctions),
    matching_function=MinimumWeightMatchingFunction(
        columns=FSMatchingColumns,
        weights=ones(7),
    ),
    minimum_area=300, # Minimum floe area for tracking
    maximum_area=90e3, # Maximum floe area for tracking
    maximum_time_step=Day(2), # Maximum length of time to skip
)
    return FloeTracker(;
        filter_function, matching_function, minimum_area, maximum_area, maximum_time_step
    )
end

end
