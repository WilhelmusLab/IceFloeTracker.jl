"""

This module documents the simplified IFT pipeline calibrated and validated in the paper
"Automatic segmentation and tracking of sea ice floes using Ice Floe Tracker algorithm" by
Watkins et al., submitted to IEEE Transactions in Geoscience and Remote Sensing (TGRS).

The pipeline contains distinct struct/functor pairs and uses both semantic and object-based 
analysis techniques. It introduces Preprocess() and Classify() functors, extending the Segment()
and Track() functors in LopezAcosta2019 and LopezAcosta2019Tiling.

"""

module TGRS2026

using Images
using DataFrames
import Dates: Day
import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering:
    ContrastLimitedAdaptiveHistogramEqualization

import ..ImageUtils:
    get_tiles

import ..Morphology:
    strel_disk

import ..Preprocessing:
    apply_landmask,
    apply_landmask!,
    Watkins2026CloudMask

import ..Segmentation:
    dist_morph_split,
    get_relevant_set,
    kmeans_segmentation,
    kmeans_binarization,
    IceDetectionBrightnessMidpoint,
    IceDetectionBrightnessPeaksMODIS721,
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

# Q: does image sharpening, nonlinear filtering change the quality of the result? nonlinear filtering 
# in particular is expensive.
function (p::Preprocess)(
    image::AbstractArray{<:Union{AbstractGray,TransparentGray,AbstractRGB,TransparentRGB}}, landmask
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

# Q: Is using both Band 1 and Band 2 necessary? 
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
    ice = Gray.(blue.(apply_landmask(fc_masked, clouds))) |> ice_mask_algorithm

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
    (max_hole_fill=1500, max_depth=10, max_depth_ratio=0.5, max_expand=3, opening_strel=strel_box((3, 3))),
    (max_hole_fill=2500, max_depth=25, max_depth_ratio=0.5, max_expand=3, opening_strel=strel_disk(3))
]
floe_filtering_params = (
    minimum_floe_size=64,
    maximum_floe_size=90e3,
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
    floe_splitting_algorithm = dist_morph_split # TODO: add binarization algorithm
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
    falsecolor_image = float64.(RGB.(falsecolor))
    coastal_buffer_mask = reinterpret(Bool, coastal_buffer_mask)
    landmask = reinterpret(Bool, landmask)

    # n, m = size(truecolor_image) TODO: warn if the tile size ends up smaller than 400 px
    tiles = get_tiles(truecolor; s.tile_settings...)

    @info "Preprocess"
    preproc_gray = s.preprocessing_algorithm(truecolor_image, landmask)

    @info "Classify"
    classifier = s.classification_algorithm
    classified_image = classifier(falsecolor_image, landmask)
    masks = Dict(k => classified_image .== classifier.key[k] for k in keys(classifier.key))
    push!(masks, "coastal_buffer_mask" => coastal_buffer_mask)

    # Then check for sufficient ocean pixels (speed up for large images)
    filtered_tiles = filter(
        t -> sum(.!masks["land"][t...]) > s.min_ocean_pixel_count, tiles
    )

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
    # final_floes = merge_floes(candidate_splits, falsecolor_image; p.floe_merging_params...)
    final_floes = candidate_splits[1]

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
        segment_mean_truecolor=n0f8.(segment_mean_map(segments_tc)) # dmw: does this work with view_seg()?
        segment_mean_falsecolor=n0f8.(segment_mean_map(segments_fc))
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            coastal_buffer_mask=Gray.(masks["coastal_buffer_mask"]),
            cloud_mask=Gray.(masks["cloud"]),
            ice_mask=Gray.(masks["ice"]),
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
function kmeans_binarization_multiclass(preproc_gray, falsecolor_image, masks, tiles;
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
        band_7_max=clear_sky_ice_params.b7,
        possible_ice_threshold=clear_sky_ice_params.b2,
        nbins=128, minimum_prominence=0.01)

    cloudy_ice_kmeans = kmeans_binarization(preproc_gray, fc_masked, tiles;
        k=cloudy_ice_params.k, cluster_selection_algorithm=cloudy_ice_detector
    ) .> 0
    clear_sky_ice_kmeans = kmeans_binarization(preproc_gray, fc_masked, tiles;
        k=clear_sky_ice_params.k, cluster_selection_algorithm=clear_sky_ice_detector
    ) .> 0

    clear_sky_ice_kmeans[masks[mask_cloud_key] .> 0] .= cloudy_ice_kmeans[masks[mask_cloud_key] .> 0]

    return clear_sky_ice_kmeans
end

"""
    extended_regionprops_table()

Calls @ref[`regionprops_table`] with the provided `properties` list. Then, adds information on
floe-average overlap with the provided `masks` (expects Dict with mask name => binary mask), 
band-average reflectance from the falsecolor image, and band 1 boundary contrast. Finally, uses
a provided probability function to add a `probability` column indicating the likelihood the object
is an ice floe.

"""
function extended_regionprops_table(
    img_indexmap,
    falsecolor_image,
    masks;
    boundary_radius=15,
    properties=[
        :label, :area, :perimeter, :bbox,
        :centroid, :convex_area, :major_axis_length,
        :minor_axis_length, :orientation,
        :circularity, :solidity],
    probability_function=LogisticRegressionFilter,
    convex_area_algorithm=PolygonConvexArea(),
)
    img_indexmap = copy(img_indexmap)
    indices = component_indices(img_indexmap)

    props_df = regionprops_table(img_indexmap;
        properties=properties,
        convex_area_algorithm=convex_area_algorithm,
    )
    # Return empty dataframe if no floes in image
    nrow(props_df) == 0 && return props_df

    transform!(props_df, :area => ByRow(x -> x^0.5) => :length_scale)

    # Don't allow circularity or solidity greater than 1
    transform!(props_df, :solidity => ByRow(x -> minimum([x, 1])) => :solidity)
    transform!(props_df, :circularity => ByRow(x -> minimum([x, 1])) => :circularity)

    # Get the average area coverage for each of the masks
    mask_mean(r, mask) = mean(mask[indices[r]])
    for k in keys(masks)
        props_df[:, Symbol(k, "_fraction")] = mask_mean.(props_df[:, :label], [masks[k]])
    end

    # Get the mean reflectance and mean boundary reflectance, and expand the results into named color channels
    add_mean_reflectance!(props_df, falsecolor_image, indices)
    add_mean_boundary_reflectance!(props_df, falsecolor_image, img_indexmap; radius=boundary_radius)

    # TODO: generalize with a map from channel number to channel name
    # Could make this a for loop with transform!()
    props_df[:, :b7_mean_reflectance] = red.(props_df.mean_reflectance)
    props_df[:, :b2_mean_reflectance] = green.(props_df.mean_reflectance)
    props_df[:, :b1_mean_reflectance] = blue.(props_df.mean_reflectance)

    props_df[:, :b7_mean_boundary_reflectance] = red.(props_df.mean_boundary_reflectance)
    props_df[:, :b2_mean_boundary_reflectance] = green.(props_df.mean_boundary_reflectance)
    props_df[:, :b1_mean_boundary_reflectance] = blue.(props_df.mean_boundary_reflectance)

    props_df[:, :b7_mean_boundary_contrast] = props_df[:, :b7_mean_reflectance] .- props_df[:, :b7_mean_boundary_reflectance]
    props_df[:, :b2_mean_boundary_contrast] = props_df[:, :b2_mean_reflectance] .- props_df[:, :b2_mean_boundary_reflectance]
    props_df[:, :b1_mean_boundary_contrast] = props_df[:, :b1_mean_reflectance] .- props_df[:, :b1_mean_boundary_reflectance]

    # TODO: generalize to include inplace option
    props_df[:, :probability] .= probability_function(props_df)

    # Drop the RGB columns in the returned dataframe
    return props_df[:, Not(:mean_reflectance, :mean_boundary_reflectance)]
end

"""LogisticRegressionFilter(df;
    coefs = Dict(
        "intercept"           => -97.1879,
        "length_scale"        => 0.1267,
        "solidity"            => 91.164,
        "b1_mean_reflectance" => 7.354,
        "b7_mean_reflectance" => -1.517,
        "b1_mean_boundary_contrast" => 2.239,
        )
    )
    LogisticRegressionFilter!(df; coefs)

Apply the logistic regression function with the provided set of coefficients. The in-place version
adds a column "probability" to the dataframe, while the non-in-place version returns a vector
with probabilities.

"""
function LogisticRegressionFilter(df;
    coefs=Dict(
        "intercept" => -97.1879,
        "length_scale" => 0.1267,
        "solidity" => 91.164,
        "b1_mean_reflectance" => 7.354,
        "b7_mean_reflectance" => -1.517,
        "b1_mean_boundary_contrast" => 2.239,
    )
)
    colnames = [x for x in keys(coefs)]
    b = [x for x in values(coefs)]
    df[:, :intercept] .= 1
    df_ = copy(df)[:, colnames]
    return 1 ./ (1 .+ exp.(-Matrix(df_[:, colnames]) * b))
end

function LogisticRegressionFilter!(df;
    coefs=Dict(
        "intercept" => -97.1879,
        "length_scale" => 0.1267,
        "solidity" => 91.164,
        "b1_mean_reflectance" => 7.354,
        "b7_mean_reflectance" => -1.517,
        "b1_mean_boundary_contrast" => 2.239,
    )
)
    colnames = [x for x in keys(coefs)]
    b = [x for x in values(coefs)]
    df[:, :intercept] = 1
    df[:, :probability] = 1 ./ (1 .+ exp.(-Matrix(df[:, colnames]) * b))
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
the bounding boxes in `props_df` so that they don't have to be re-computed.
"""
function add_mean_boundary_reflectance!(props_df, img, labels; radius=15)
    n, m = size(labels)
    bdry_ref = []
    for data in eachrow(props_df)
        # expand the bounding box by radius
        # minimum row is the maximum 
        rmin = maximum((data.min_row - radius, 1))
        rmax = minimum((data.max_row + radius, n))
        cmin = maximum((data.min_col - radius, 1))
        cmax = minimum((data.max_col + radius, m))

        label_subset = Int64.(labels[rmin:rmax, cmin:cmax] .== data.label)
        boundary = expand_labels(label_subset, radius)
        boundary[label_subset .> 0] .= 0
        image_subset = img[rmin:rmax, cmin:cmax]
        push!(bdry_ref, mean(image_subset[boundary .> 0]))
    end
    props_df.mean_boundary_reflectance = bdry_ref
end


"""
    compare_objects(
        df1, df2, labels1, labels2;
        indices1=component_indices(labels1),
        indices2=component_indices(labels2),
        comp_properties=[
            :label, :area, :row_centroid, :col_centroid,
            :max_col, :max_row, :min_col, :min_row, :probability
        ],
        tol_area_fraction=0.05,
    )

Produce a dataframe comparing objects in a pair of labeled images, including 
all paired labels between labels1 and labels2 with area overlap greater than
`tol_area_fraction` relative to either label. Additionally computes the distance between centroids, area overlap, and fractional area overlap.

Inputs:
    - `df1` = region properties dataframe from labels1
    - `df2` = region properties dataframe from labels2
    - `labels1` = labeled image (Matrix{Int64})
    - `labels2` = labeled image (Matrix{Int64})
    - `indices1=component_indices(labels1)` = Indices map, option to reuse from earlier in processing 
    - `indices2=component_indices(labels2)` = Indices map, option to reuse from earlier
    - `comp_properties=[
            :label, :area, :row_centroid, :col_centroid,
            :max_col, :max_row, :min_col, :min_row, :probability
        ]` = Columns in df1 and df2 to include in comparison
    - `tol_area_fraction=0.05`= Minimum area fraction to include in comparison
"""
function compare_objects(
    df1::DataFrame,
    df2::DataFrame,
    labels1::Matrix{Int64},
    labels2::Matrix{Int64}; # Should this be keyword or no?
    indices1=component_indices(labels1),
    indices2=component_indices(labels2),
    comp_properties=[
        :label, :area, :row_centroid, :col_centroid,
        :max_col, :max_row, :min_col, :min_row, :probability
    ],
    tol_area_fraction=0.05, # TODO: decide whether we should filter probability here
)::DataFrame

    # Get list of labels in 1 with nonzero intersection
    no_overlaps1 = _nonoverlapping_labels(labels2, indices1, df1.label)
    overlaps1 = setdiff(df1.label, no_overlaps1)

    # Make list of intersections from 1 to 2
    s1_label_list = []
    s2_label_list = []
    for r in overlaps1
        for s in filter(r -> r != 0, unique(labels2[indices1[r]]))
            append!(s1_label_list, r)
            append!(s2_label_list, s)
        end
    end

    # Generate joint dataframe
    df_comp1 = rename(df1[:, comp_properties],
        Dict(p => Symbol("s1_", p) for p in comp_properties))
    df_comp2 = rename(df2[:, comp_properties],
        Dict(p => Symbol("s2_", p) for p in comp_properties))
    df_dict1 = Dict(row.s1_label => row for row in eachrow(df_comp1))
    df_dict2 = Dict(row.s2_label => row for row in eachrow(df_comp2))
    df_comp = hcat(
        DataFrame([df_dict1[l] for l in s1_label_list]),
        DataFrame([df_dict2[l] for l in s2_label_list])
    )

    # Compute overlap metrics
    transform!(df_comp,
        [:s1_row_centroid, :s2_row_centroid,
            :s1_col_centroid, :s2_col_centroid] =>
            ByRow((r1, r2, c1, c2) -> sqrt((r1 - r2)^2 + (c1 - c2)^2)) =>
                :s1_s2_dist
    )

    transform!(df_comp,
        [:s1_label, :s2_label,
            :s1_min_row, :s1_max_row, :s1_min_col, :s2_max_col] =>
            ByRow((l1, l2, rmin, rmax, cmin, cmax) ->
                sum(
                    (labels1[rmin:rmax, cmin:cmax] .== l1) .&&
                        (labels2[rmin:rmax, cmin:cmax] .== l2)
                )
            ) =>
                :s1_s2_area_overlap
    )

    transform!(df_comp,
        [:s1_s2_area_overlap, :s1_area] => ByRow((a0, a1) -> a0/a1) =>
            :s1_area_fraction
    )

    transform!(df_comp,
        [:s1_s2_area_overlap, :s2_area] => ByRow((a0, a1) -> a0/a1) =>
            :s2_area_fraction
    )

    subset!(df_comp, :s1_area_fraction => r -> r .> tol_area_fraction)
    subset!(df_comp, :s2_area_fraction => r -> r .> tol_area_fraction)

    return df_comp
end

"""
    _nonoverlapping_labels(other, indices, labels)

Return a list of labels in matrix `other` which have no
overlap with the list of labels `labels` and the corresponding
indices dictionary `indices`. Both `labels` and `indices` 
come from a second labeled indexmap to be compared with `other`.

"""
function _nonoverlapping_labels(other, indices, labels)
    return [
        label for label in labels
                  if maximum(other[indices[label]]) == 0
    ]
end

"""
    _assign_labels!(output, indices, labels; offset=0)

Insert each label from list `labels` into `output` using 
the indices dictionary `indices`. Optional `offset` integer
can be added to avoid duplicating an existing label.

"""
function _assign_labels!(output, indices, labels; offset=0)
    foreach(labels) do label
        output[indices[label]] .= label + offset
    end
end

"""
    _remove_labels!(output, indices, remove_labels)

Remove regions of `output` by setting the indices to 0.
The labels in `remove_labels` correspond to the dictionary
keys in `indices`.

"""
function _remove_labels!(output, indices, remove_labels)
    foreach(remove_labels) do label
        output[indices[label]] .= 0
    end
end

"""
    merge_arrays!(
        output,
        indices1,
        indices2,
        remove_labels,
        add_labels
    )

Update output by (1) removing the labels for each `L1` the list `remove_labels` by setting everything in `indices1[L1]` to 0 and then (2)
writing `L2` into labels1 for each `L2` in `add_labels`.
"""
function merge_arrays!(output, indices1, indices2, remove_labels, add_labels)
    _remove_labels!(output, indices1, remove_labels)
    _assign_labels!(output, indices2, add_labels;
        offset=maximum(labels1))
end


"""
    colorize_classification(labeled_image; color_map)

Convenience function to make an image from a classified image, using
the dictionary `color_map=Dict(label_integer => color)`.
"""
function colorize_classification(labeled_image;
    color_map=Dict(
        0=>RGB(0),
        1=>RGB(0.018, 0.49, 0.64),
        2=>RGB(1),
        3=>RGB(0.84, 0.73, 0.94)
    )
)
    return n0f8.(map(i -> color_map[i], labeled_image))
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

const FilterFunctions = [
    max_travel_distance_filter,
    area_relative_error_filter,
    convex_area_relative_error_filter,
    major_axis_relative_error_filter,
    minor_axis_relative_error_filter,
    shape_difference_filter,
    psi_s_correlation_filter,
]

const MatchingColumns = [
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
function Track(;
    filter_function=ChainedFilterFunction(; filters=FilterFunctions),
    matching_function=MinimumWeightMatchingFunction(
        columns=MatchingColumns,
        weights=ones(7),
    ),
    minimum_area=300, # Minimum floe area for tracking
    maximum_area=90e3, # Maximum floe area for tracking
    maximum_time_step=Day(1), # Maximum length of time to skip
)
    return FloeTracker(;
        filter_function, matching_function, minimum_area, maximum_area, maximum_time_step
    )
end

end
