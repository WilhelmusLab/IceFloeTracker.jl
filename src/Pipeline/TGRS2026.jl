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

#### Segmentation parameters
# Potential changes: add a Binarize() method to be applied to the preprocessed image.
tile_settings = (; rblocks=1, cblocks=1)
min_ocean_pixel_count = 5000
preprocessing_algorithm = Preprocess()
classification_algorithm = Classify()
floe_splitting_params = [
        (max_hole_fill=500, max_depth=5, max_depth_ratio=0.5, max_expand=3, opening_strel=strel_disk(1)),
        (max_hole_fill=1500, max_depth=10, max_depth_ratio=0.5, max_expand=3, opening_strel=strel_box((3,3))),
        (max_hole_fill=2500, max_depth=25, max_depth_ratio=0.5, max_expand=3, opening_strel=strel_disk(3))
    ]
floe_filtering_params = (
    minimum_floe_size = 64,
    maximum_floe_size = 90e3,
    )
"""
    Segment()

Function to identify pixel types and detect ice floes from MODIS true color and false color imagery.
The coastal buffer mask is used to identify potential landfast ice segments. 

- `coastal_buffer_structuring_element`: Structuring element to dilate the coast mask. Used to identify landfast ice.
- `tile_settings`: Named tuple with empty first argument to provide keywords. Default is (; rblocks=1, cblocks=1)  
- `min_ocean_pixel_count`: Minimum number of non-land pixels to run the algorithm on a tile.
- `preprocessing_algorithm`: Algorithm to generate the grayscale image to send to the binarization stage.
- `classification_algorithm`: Algorithm to identify land, ice, water, and clouds
- `floe_splitting_algorithm`: Currently only tested with the dist morph split algorithm. The function will apply the 
    splitting algorithm for each set of parameters supplied in the next argument.
- `floe_splitting_params`: List of named tuples with parameters for the `floe_splitting_algorithm`.
- `floe_filtering_algorithm`: Algorithm to remove likely non-floes from the merged segmentation result.

"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_disk(25)
    tile_settings = tile_settings
    min_ocean_pixel_count = min_ocean_pixel_count
    preprocessing_algorithm = preprocessing_algorithm
    classification_algorithm = classification_algorithm
    floe_splitting_algorithm = dist_morph_split
    floe_splitting_params = floe_splitting_params
    floe_filtering_params = filter_floes_params
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
    coastal_buffer_mask = reinterpret(Bool, coastal_buffer_mask)
    landmask = reinterpret(Bool, landmask)

    # n, m = size(truecolor_image) TODO: warn if the tile size ends up smaller than 400 px
    tiles = get_tiles(truecolor; p.tile_settings...)
    
    @info "Preprocess"
    preproc_gray = p.preprocessing_algorithm(truecolor_image, landmask)

    @info "Classify"
    classifier = p.classification_algorithm
    classified_image = classifier(falsecolor_image, landmask)
    masks = Dict(k => classified_image .== classifier.key[k] for k in keys(classifier))
    push!(masks, "coastal_buffer_mask" => coastal_buffer_mask)

    # Then check for sufficient ocean pixels (speed up for large images)
    filtered_tiles = filter(
        t -> sum(.!masks["land"][t...]) > s.min_ocean_pixel_count, tiles
    );

    @info "Detect Floes"
    binarized_image = kmeans_binarization_multiclass(
        preproc_gray, falsecolor_image, masks, filtered_tiles,
    )
    candidate_splits = [
        dist_morph_split(binarized_image; pset...) for pset in s.floe_splitting_params
        ]
    
    # Size-based filter
    remove_small_segments!.(candidate_splits, s.floe_filtering_params.minimum_floe_size)
    remove_large_segments!.(candidate_splits, s.floe_filtering_params.maximum_floe_size)
    
    
    @info "Joining segmentation results"
    final_floes = merge_floes(kmeans_split_floes, adaptive_split_floes, preproc_gray)

    # Repeat size-based filter in case artifacts were created
    remove_small_segments!(final_floes, s.floe_filtering_params.minimum_floe_size)
    remove_large_segments!(final_floes, s.floe_filtering_params.maximum_floe_size)

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

#### Helper functions for segmentation ####
"""
    kmeans_binarization_multiclass(preproc_gray, falsecolor_image, masks;
    cloudy_ice_params=(k=3, b7=0.7, b2=0.56),
    clear_sky_ice_params=(k=4, b7=0.18, b2=0.46),
    mask_land_key="coastal_buffer_mask",
    mask_cloud_key="cloud",
    )

Wrapper for the kmeans binarization function. Uses one set of k-means settings for clear
sky scenes and another set for cloud-covered scenes. The parameters are used to initialize
the IceDetectionBrightnessPeaksMODIS721 algorithm.
"""
function kmeans_binarization_multiclass(preproc_gray, falsecolor_image, masks;
    cloudy_ice_params=(k=3, b7=0.7, b2=0.56),
    clear_sky_ice_params=(k=4, b7=0.18, b2=0.46),
    mask_land_key="coastal_buffer_mask",
    mask_cloud_key="cloud",
    )

    fc_masked = apply_landmask(falsecolor_image, masks[mask_land_key])
    
    cloudy_ice_detector = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=cloudy_ice_params.b7,
        possible_ice_threshold=cloudy_ice_params.b2,
        nbins=128, minimum_prominence=0.03)
    clear_sky_ice_detector = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=clearsky_ice_params.b7,
        possible_ice_threshold=clearsky_ice_params.b2,
        nbins=128, minimum_prominence=0.01);

    cloudy_ice_kmeans = kmeans_binarization(preproc_gray, fc_masked;
        k=cloudy_ice_params.k, cluster_selection_algorithm=cloudy_ice_detector
    ) .> 0
    clear_sky_ice_kmeans = kmeans_binarization(preproc_gray, fc_masked;
        k=clear_sky_ice_params.k, cluster_selection_algorithm=clear_sky_ice_detector
    ) .> 0
    
    clear_sky_ice_kmeans[masks[mask_cloud_key] .> 0] .= cloudy_ice_kmeans[masks[mask_cloud_key] .> 0]

    return clear_sky_ice_kmeans
end


function extended_regionprops(
    img_indexmap,
    coastal_buffer_mask,
    classified_image,
    falsecolor_image; # expects band 7-2-1
    boundary_radius=15,
    classification_key=Dict("land"=>0, "water"=>1, "ice"=>2, "cloud"=>3),
    properties = [
        :label, :area, :perimeter, :bbox,
        :centroid, :convex_area, :major_axis_length,
        :minor_axis_length, :orientation,
        :circularity, :solidity],
    probability_function=LogisticRegressionFilter,
)
    img_indexmap = copy(img_indexmap)
    indices = component_indices(img_indexmap)

    results_df = regionprops_table(img_indexmap;
        properties=properties,
        convex_area_algorithm=PolygonConvexArea()
    )
    # Return blank image if no floes remain
    nrow(results_df) == 0 && return results_df

    results_df[:, :length_scale] = results_df[:, :area] .^ 0.5
    # Correct circularity error
    results_df[:, :circularity] = 4 * pi * results_df[:, :area] ./ results_df[:, :perimeter] .^ 2
    
    mask_mean(r, mask) = mean(mask[indices[r]])
    masks = Dict(k => classified_image .== classification_key[k] for k in keys(classification_key))
    push!(masks, "coast" => coastal_buffer_mask)
    
    results_df[:, :cloud_fraction] =  mask_mean.(results_df[:, :label], [masks["cloud"]])
    results_df[:, :ice_fraction] =  mask_mean.(results_df[:, :label], [masks["ice"]])
    results_df[:, :water_fraction] =  mask_mean.(results_df[:, :label], [masks["water"]])
    results_df[:, :coastal_buffer_fraction] =  mask_mean.(results_df[:, :label], [masks["coast"]])
    
    # Compute mean reflectance, assuming the input image is MODIS False Color
    segment_mean_reflectance = Dict(r => mean(falsecolor_image[indices[r]]) for r in keys(indices))
    b = (r -> segment_mean_reflectance[r]).(results_df[:, :label])
    results_df[:, :b1_reflectance_mean] = blue.(b)
    results_df[:, :b7_reflectance_mean] = red.(b)
    results_df[:, :b2_reflectance_mean] = green.(b)

    # Compute mean Band 1 boundary reflectance using the boundary radius for an expansion limit
    b1 = blue.(falsecolor_image)
    eroded_labels = img_indexmap .* erode(img_indexmap .> 0)
    bdry_indexmap = expand_labels(img_indexmap, boundary_radius) .- eroded_labels
    bdry_indices = component_indices(bdry_indexmap)
    bdry_labels = intersect(results_df[:, :label], unique(bdry_indexmap))
    b1_bdry_means = Dict(L => mean(b1[bdry_indices[L]]) for L in bdry_labels)
    for L ∈ results_df[:, :label]
        if L ∉ bdry_labels
            push!(b1_bdry_means, L => 0)
        end
    end
    results_df[:, :b1_reflectance_bdry_mean] = [b1_bdry_means[L] for L in results_df[:, :label]]
    results_df[:, :b1_bdry_contrast] = results_df[:, :b1_reflectance_mean] .- results_df[:, :b1_reflectance_bdry_mean]
    
    results_df[:, :probability] .= probability_function(results_df)
    return results_df
end

"""LogisticRegressionFilter(df;
    coefs = Dict(
        "intercept"           => -97.1879,
        "length_scale"        => 0.1267,
        "solidity"            => 91.164,
        "b1_reflectance_mean" => 7.354,
        "b1_bdry_contrast"    => 2.239,
        "b7_reflectance_mean" => -1.517,
        )
    )
    LogisticRegressionFilter!(df; coefs)

Apply the logistic regression function with the provided set of coefficients. The in-place version
adds a column "probability" to the dataframe, while the non-in-place version returns a vector
with probabilities.

"""
function LogisticRegressionFilter(df;
    coefs = Dict(
        "intercept"           => -97.1879,
        "length_scale"        => 0.1267,
        "solidity"            => 91.164,
        "b1_reflectance_mean" => 7.354,
        "b1_bdry_contrast"    => 2.239,
        "b7_reflectance_mean" => -1.517,
        )
    )
    colnames = [x for x in keys(coefs)]
    b = [x for x in values(coefs)]
    df_ = copy(df)[:, colnames]
    df_[:, :intercept] .= 1;
    return 1 ./ (1 .+ exp.(-Matrix(df_[:, colnames]) * b))
end

function LogisticRegressionFilter!(df;
    coefs = Dict(
        "intercept"           => -97.1879,
        "length_scale"        => 0.1267,
        "solidity"            => 91.164,
        "b1_reflectance_mean" => 7.354,
        "b1_bdry_contrast"    => 2.239,
        "b7_reflectance_mean" => -1.517,
        )
    )
    colnames = [x for x in keys(coefs)]
    b = [x for x in values(coefs)]
    df[:, :intercept] .= 1;
    df[:, :probability] = 1 ./ (1 .+ exp.(-Matrix(df_[:, colnames]) * b))
end

"""
    add_mean_reflectance!(props_df, img, indices)

Compute the mean reflectance for `img` for each label in `props_df`. Assumes
that `props_df` contains labels corresponding to the dictionary `indices` 
(see @ref[`component_indices`]).
"""
function add_mean_reflectance!(props_df, img, indices)
    segment_mean_reflectance(r) = mean(img[indices[r]])
    props_df.mean_reflectance = segment_mean_reflectance.(props_df.label)
end

"""
    add_mean_boundary_reflectance!(props_df, img, labels; radius=15)

Compute the average of `img` within `radius` of the objects in `labels`. Uses
the bounding boxes in `regionprops_df` so that they don't have to be re-computed.
"""
function add_mean_boundary_reflectance!(props_df, img, labels; radius=15)
    n, m = size(labels)
    bdry_ref = zeros(Float64, nrow(regionprops_df))
    for data in eachrow(regionprops_df)
        # expand the bounding box by radius
        # minimum row is the maximum 
        rmin = maximum((data.min_row - radius, 0))
        rmax = mimimum((data.max_row + radius, n))
        cmin = maximum((data.min_col - radius, 0))
        cmax = minimum((data.max_col + radius, m))

        label_subset = Int64.(labels[rmin:rmax, cmin:cmax] .== data.label)
        boundary = expand_labels(label_subset, radius)
        boundary[label_subset .> 0] .= 0
        image_subset = img[rmin:rmax, cmin:cmax]
        push!(bdry_ref, mean(vec(image_subset[boundary .> 0])))
    end
    props_df.mean_boundary_reflectance = bdry_ref
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
