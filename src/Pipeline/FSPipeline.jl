module FSPipeline
"""

Simplified segmentation pipeline with calibrated parameters for the Greenland Sea / Fram Strait workflow.

"""

using Images
using DataFrames
import Dates: Day
import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering:
    nonlinear_diffusion,
    PeronaMalikDiffusion,
    unsharp_mask,
    ContrastLimitedAdaptiveHistogramEqualization

import ..Morphology: fill_holes, strel_disk

import ..Preprocessing:
    create_landmask,
    create_cloudmask,
    apply_landmask,
    apply_landmask!,
    apply_cloudmask,
    apply_cloudmask!,
    Watkins2026CloudMask

import ..ImageUtils: get_tiles, imbrighten
import ..Segmentation:
    component_perimeters,
    expand_labels,
    get_relevant_set,
    kmeans_binarization,
    tiled_adaptive_binarization,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionBrightnessMidpoint,
    PolygonConvexArea,
    regionprops_table,
    remove_small_segments!,
    remove_large_segments!,
    segment_mean_map,
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

# Preprocess Params
diffusion_algorithm = PeronaMalikDiffusion(; λ=0.1, K=0.1, niters=7, g="exponential")
adapthisteq_params = (nbins=256, rblocks=4, cblocks=4, clip=1)
unsharp_mask_params = (radius=50, amount=0.3, threshold=0.01)

"""
   Preprocess(
        diffusion_algorithm = PeronaMalikDiffusion(λ=0.1, K=0.1, niters=5, g="exponential")
        adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.99) # rblocks/cblocks not used yet -- add with CLAHE.jl
        unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
    )
    Preprocess()(img, cloudmask, landmask)

    Converts input image to grayscale, then preprocesses by appling nonlinear diffusion,
    adaptive histogram equalization, and unsharp masking. Diffusion and unsharp masking are applied
    to each tile, while the adaptive histogram equalization is divided according to the parameter
    specifications.

    Note: results are strongly sensitive to the choice of rblocks, cblocks, and clipping. Large clipping parameters with
    small blocks results in noisy images and poor performance. With larger blocks, a higher clipping parameter can help.

"""
@kwdef struct Preprocess <: IceFloePreprocessingAlgorithm
    diffusion_algorithm = diffusion_algorithm
    adapthisteq_params = adapthisteq_params
    unsharp_mask_params = unsharp_mask_params
end

function (p::Preprocess)(
    truecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}}, landmask, tiles
)
    # Cast to grayscale first to save compute time
    proc_img = Gray.(truecolor_image)

    # Diffusion and sharpening
    proc_img .= nonlinear_diffusion(proc_img, tiles, p.diffusion_algorithm)

    adjust_histogram!(
        proc_img,
        ContrastLimitedAdaptiveHistogramEqualization(;
            nbins=p.adapthisteq_params.nbins,
            rblocks=p.adapthisteq_params.rblocks,
            cblocks=p.adapthisteq_params.cblocks,
            clip=p.adapthisteq_params.clip,
        ),
    )

    proc_img .= unsharp_mask(
        proc_img,
        p.unsharp_mask_params.radius,
        p.unsharp_mask_params.amount,
        p.unsharp_mask_params.threshold,
    )

    # Re-apply mask so sharpening doesn't bleed into land
    apply_landmask!(proc_img, landmask)
    return proc_img
end

# Default segmentation parameters
coastal_buffer_structuring_element = strel_box((51, 51))
cloud_mask_algorithm = Watkins2026CloudMask()
preprocessing_algorithm = Preprocess()
tile_size_pixels = 400
min_tile_ice_pixel_count=300
preliminary_ice_mask = IceDetectionBrightnessMidpoint(; minimum_reflectance=0.3)
kmeans_params = (
    k=4,
    maxiter=50,
    random_seed=45,
    cluster_selection_algorithm=IceDetectionBrightnessPeaksMODIS721(;
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01,
    ),
)
adaptive_params = (window_size=400, percentage=0)
cleanup_binary_params = (
    erosion_strel=strel_box((3, 3)), init_max_fill=100, conditional_max_fill=500
)
floe_splitting_params = (max_hole_fill=2000, max_distance=5, max_expand=3)
floe_filtering_params = (
    min_floe_size=100,
    max_floe_size=90_000,
    boundary_radius=15,
    min_circularity=0.3,
    min_reflectance=0.4,
    min_contrast=0.01,
    min_probability=0.5,
)
floe_merging_params = (
    distance_threshold_pixels=10, area_error_threshold=0.25, min_floe_size=100
)

"""
    FSPipeline.Segment()

Segmentation routine for identifying moderate to large floes in the Fram Strait.
The image preprocessing is supplied as an function in the functor setup.


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
    tile_size_pixels = tile_size_pixels
    min_tile_ice_pixel_count = min_tile_ice_pixel_count
    preliminary_ice_mask = preliminary_ice_mask
    kmeans_params = kmeans_params
    adaptive_params = adaptive_params
    cleanup_binary_params = cleanup_binary_params
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
    land_mask = landmask .> 0 # make sure it's a bitmatrix, in case it's passed as Gray
    apply_landmask!(truecolor_image, land_mask)
    apply_landmask!(falsecolor_image, land_mask)

    n, m = size(truecolor_image)
    tile_size_pixels = s.tile_size_pixels
    nmin, nmax = extrema([n, m])
    tile_size_pixels > nmax && begin
        @warn "Tile size too large; clamping to min(height, width)."
        tile_size_pixels = nmin
    end

    (nr, nc) = round.(Int, size(truecolor_image) ./ tile_size_pixels)
    tiles = get_tiles(truecolor_image; rblocks=nr, cblocks=nc)

    @info "Building masks"
    cloud_mask = create_cloudmask(falsecolor_image, s.cloud_mask_algorithm)

    # 2. Intermediate images - apply coastal buffer and cloud mask
    joint_mask = coastal_buffer_mask .|| cloud_mask
    tc_masked = apply_landmask(truecolor_image, joint_mask)
    fc_masked = apply_landmask(falsecolor_image, joint_mask)

    # First check for sufficient non-land and non-cloud pixels
    filtered_tiles = filter(
        t -> sum(.!joint_mask[t...]) > s.min_tile_ice_pixel_count, tiles
    );

    # Then check for sufficient possible sea ice pixels
    prelim_ice_mask = s.preliminary_ice_mask(Gray.(red.(tc_masked)), filtered_tiles)
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > s.min_tile_ice_pixel_count, filtered_tiles
    );
    water_mask = .!(prelim_ice_mask .|| cloud_mask .|| land_mask)

    @info "Preprocessing truecolor image"
    preproc_gray = float64.(
        s.preprocessing_algorithm(truecolor_image, land_mask, filtered_tiles)
    );

    @info "Segmentation"

    adaptive_result = apply_landmask(
        binarize(preproc_gray, AdaptiveThreshold(; s.adaptive_params...)) .> 0,
        water_mask .|| land_mask
    )

    kmeans_result = kmeans_binarization(
        preproc_gray, fc_masked, filtered_tiles; s.kmeans_params...
    )

    clean_and_split(r) = dist_morph_split(
            clean_binary_floes(
                r, prelim_ice_mask, cloud_mask; s.cleanup_binary_params...
            );
        s.floe_splitting_params...,
    )

    labeled_images = clean_and_split.([kmeans_result, adaptive_result])
    
    @info "Filter and merge"
    # Update the filtering method so that the regionprops table is only called once. Now, it returns 
    # a dataframe so that the merge floes can take two tuples as inputs.
    filtered_floes = filter_floes.(labeled_images, [coastal_buffer_mask], [cloud_mask], [falsecolor_image]; s.floe_filtering_params...)
    keep_labels!.(labeled_images, filtered_floes .|> r -> r.label)   

    # TODO: Update merge floes to just use the properties in the filter floes table
    final_floes = merge_floes(filtered_floes..., labeled_images...)

    # Remove any stray segments left over from the merge function
    remove_small_segments!(final_floes, s.floe_filtering_params.min_floe_size)
    
    # Re-label so there are no missing numbers in the component list

    final_floes .= label_components(final_floes)

    # Generate images with object-average colors
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
            kmeans_binarized=kmeans_result .> 0,
            adaptive_binarized=adaptive_result .> 0,
            kmeans_floes=labeled_images[1] .> 0,
            adaptive_floes=labeled_images[2] .> 0,
            final_floes=colorview_random,
            labels_map=final_floes,
            segment_mean_falsecolor=segment_mean_falsecolor,
            segment_mean_truecolor=segment_mean_truecolor,
        )
    end
    return segments_tc
end

"""
    clean_binary_floes(binary_img, icemask, cloudmask;
        erosion_strel=strel_box((3,3)),
        init_max_fill=100,
        conditional_max_fill=500
    )

Fill holes in a binary mask. First, fill holes in the eroded floe shapes
up to size `init_max_fill`. Then, fill holes up to `conditional_max_fill`
if those holes are either ice or cloud. Finally, reset any filled holes that
intersect with the boundary.

"""
function clean_binary_floes(
    binary_img,
    icemask,
    cloudmask;
    erosion_strel=strel_box((3, 3)),
    init_max_fill=100,
    conditional_max_fill=500,
)
    out = deepcopy(binary_img)
    # 1. Shrink objects using the provided structuring element
    eroded_img = erode(out, erosion_strel)

    # 2. After shrinking, fill holes
    filled = .!imfill(.!eroded_img, (0, init_max_fill)) # Test how permissive this is. Should we use imfill instead?

    # 3. Identify filled holes which are part of the ice mask or the cloud mask
    filled .= filled .&& (icemask .|| cloudmask)
    filled .= .!imfill(.!filled, (0, conditional_max_fill))

    # 4. Use morphological closing to further limit openings
    filled .= closing(filled, erosion_strel)

    # 5. Set any of these filled pixels to 1 in the output image.
    out[filled .> 0] .= 1
    opening!(out)

    # 6. If the filled region intersects with a boundary, remove it
    filled .= out .!= binary_img
    out[filled .&& .! clearborder(filled)] .= 0

    return out
end

"""
    dist_morph_split(
        binary_floes::BitMatrix;
        min_floe_size::Int64=64,
        max_hole_fill::Int64=2000,
        max_distance::Int64=5,
        max_expand::Int64=3,
        strel=strel_disk(3)
    )

Method to split objects in a binary image using image morphology and the distance transform. The algorithm
operates by calculating the distance transform, which computes the distance from each labeled pixel to the background.
There are two steps: creating a ``pyramid'', then stepping down from the top of the pyramid and re-labeling or expanding
shapes as needed.

For each distance d up to `max_distance`, select pixels that are greater than that distance. Perform morphological opening,
fill holes up to `max_hole_fill`, then label components. Each of these layers is a level in the pyramid.

Then, starting from the highest level of the pyramid, check to see whether objects in the next layer down contain multiple
objects in the current layer. If an object at layer ``d-1`` contains only object at layer ``d``, then keep the object at layer ``d-1``.
Otherwise, expand the labels by `max_expand`, then intersect the expanded labels with the containing object at layer ``d-1``.

After traversing the pyramid, relabel matrix, and remove any objects smaller than the `min_floe_size`.

"""
function dist_morph_split(
    binary_floes::BitMatrix;
    max_hole_fill::Int64=2000,
    max_distance::Int64=5,
    max_expand::Int64=3,
    opening_strel=strel_disk(3),
)
    dist = distance_transform(feature_transform(.!binary_floes))
    levels = Dict(0 => label_components(opening(dist .> 0, opening_strel))) # Initialize with one run of opening
    ### Build pyramid - each size is the opened and filled thresholded image
    for dist_threshold in 0:max_distance
        markers = opening(dist .> dist_threshold, opening_strel)
        markers .= .!imfill(.!markers, (0, max_hole_fill))
        levels[dist_threshold] = label_components(markers)
    end
    final_labels = deepcopy(levels[max_distance])

    ### Descend pyramid
    for dist_threshold in max_distance:-1:1
        # Get indices from level d-1
        indices = component_indices(levels[dist_threshold - 1])

        # Expand indices at level d
        expanded = expand_labels(levels[dist_threshold], max_expand)
        for L in keys(indices)
            (L <= 0) && continue

            matched_labels = unique(levels[dist_threshold][indices[L]])

            # If intersection of the label at level
            if (0 ∈ matched_labels) && (length(matched_labels) <= 2)
                final_labels[indices[L]] .= L
                continue
            end
            # Otherwise, expand the current level, and set the next level down to the expanded indices.
            # May need to check the number of matched labels in the expanded image.
            levels[dist_threshold - 1][indices[L]] .= expanded[indices[L]]
            final_labels[indices[L]] .= expanded[indices[L]]
        end
    end
    return label_components(final_labels)
end

# Helper function for creating a filtered version of the image indexmap
# TODO: Unify approach with remove_small_segments
"""
    keep_labels!(img_indexmap, labels_list)

Given an image indexmap `img_indexmap` and a list of labels `labels_list`, 
remove any segments not in the list. In place version of `keep_labels`.

"""
function keep_labels!(img_indexmap, labels_list)    
    indices = component_indices(img_indexmap)
    labels = filter(r -> r > 0, unique(img_indexmap))
    for L in labels
        if L ∉ labels_list
            img_indexmap[indices[L]] .= 0
        end
    end
end

# TODO: dmw -- This function can be made much shorter by adding segment mean color, circularity, potentially boundary contrast to the region props function 
"""
    filter_floes(
        img_indexmap,
        coastal_buffer_mask,
        cloud_mask,
        falsecolor_image;
        min_floe_size=100,
        max_floe_size=90_000,
        expand_radius=15,
        filter_function=LogisticFilterFunction # needs to operate on a dataframe
)

Filter the image indexmap using object-wise properties. Removes objects which overlap the coastal buffer mask,
exceed the size limits, and have too-low circularity. Then, we apply the filter function to the dataframe
which by default uses a pre-fitted logistic regression function. Returns a dataframe with floe properties.

"""
function filter_floes(
    img_indexmap,
    coastal_buffer_mask,
    cloud_mask,
    falsecolor_image;
    min_floe_size=100,
    max_floe_size=90_000,
    boundary_radius=15,
    min_reflectance=0.4,
    min_circularity=0.3,
    min_solidity=0.7,
    min_contrast=0.01,
    filter_function=LogisticRegressionFilter,
    min_probability=0.5,
)
    # 1. Remove objects which overlap the coastal mask
    overlap = unique(img_indexmap[coastal_buffer_mask])
    indices = component_indices(img_indexmap)
    for L in overlap
        img_indexmap[indices[L]] .= 0
    end

    # 2. Remove objects outside the specified size bounds prior to extracting features.
    # This is important since the small features can cause problems in some feature
    # descriptors.
    remove_small_segments!(img_indexmap, min_floe_size)
    remove_large_segments!(img_indexmap, max_floe_size)

    # 3. Get object-wise properties
    results_df = regionprops_table(img_indexmap;
        properties=[:label, :area, :perimeter, :bbox, :centroid, :convex_area,
                    :major_axis_length, :minor_axis_length, :orientation],
        convex_area_algorithm=PolygonConvexArea()
    )
    # Return blank image if no floes remain
    nrow(results_df) == 0 && return results_df

    results_df[:, :length_scale] = results_df[:, :area] .^ 0.5
    results_df[:, :circularity] = 4 * π * results_df[:, :area] ./ results_df[:, :perimeter] .^ 2
    subset!(results_df, :circularity => r -> r .> min_circularity)
    results_df[:, :solidity] = results_df[:, :area] ./ results_df[:, :convex_area]
    subset!(results_df, :solidity => r -> r .> min_solidity)
    nrow(results_df) == 0 && return results_df

    results_df[:, :cloud_fraction] =  (r -> mean(cloud_mask[indices[r]])).(results_df[:, :label])
    
    # mean reflectance
    segment_mean_reflectance = segment_mean(SegmentedImage(falsecolor_image, img_indexmap))
    b = [segment_mean_reflectance[L] for L in  results_df[:, :label]]
    results_df[:, :b1_reflectance_mean] = blue.(b)
    results_df[:, :b7_reflectance_mean] = red.(b)
    results_df[:, :b2_reflectance_mean] = green.(b)
    
    subset!(results_df, :b1_reflectance_mean => r -> r .> min_reflectance)
    nrow(results_df) == 0 && return results_df
    
    # mean boundary reflectance
    b1 = blue.(falsecolor_image)
    bdry_indexmap = expand_labels(img_indexmap, boundary_radius) .- img_indexmap
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
    subset!(results_df, :b1_bdry_contrast => r -> r .> min_contrast)
    nrow(results_df) == 0 && return results_df

    results_df[:, :probability] .= filter_function(results_df)
    subset!(results_df, :probability => r -> r .> min_probability)

    return results_df
end

# TODO: dmw -- this could be a struct / functor pair, so the filter takes a FloeProbabilityFunction rather than needing to be LogisticRegression
"""
    LogisticRegressionFilter(df; coefs)

Compute the probability for each DataFrameRow using the logistic function
with coefficients defined in `coefs`. Names should include "intercept" and
names of columns in `df`.

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
    df[:, :intercept] .= 1;
    return 1 ./ (1 .+ exp.(-Matrix(df[:, colnames]) * b))
end

const default_properties = [:label, :area, :perimeter, :centroid]

"""
    objectwise_compare_segmentation(df1, df2, labels1, labels2)

Compare segmentation results using region properties. Identifies the relevant set of
objects in labels2 for each object in labels1. Expects df1 and df2 to originate
from the `filter_floes` function, so they should include columns for area, centroid,
and bounding box already. Returns a dataframe with a row for each comparison between 
labels1 and labels2 and adds comparison metrics for the distance between centroids `dist_s1_s2`
and the ratio of absolute area difference to summed area `scaled_relative_error_area_s1_s2`.

"""
function objectwise_compare_segmentation(
    df1, df2, labels1, labels2
)    
    properties = union(propertynames(df1), propertynames(df2))
    relevant_set = get_relevant_set(df1, df2, labels1, labels2)
    results = DataFrame[]
    for floe in eachrow(df1)
        g = floe.label
        g in keys(relevant_set) && begin
            df_rs = subset(df2, :label => ByRow(s -> s in relevant_set[g]))
            df_rs[:, :dist_s1_s2] = euclidean_distance(floe, df_rs; r=1) # r=1 means use pixel units, not meters
            df_rs[:, :scaled_relative_error_area] =
                abs.(df_rs.area .- floe.area) ./ (df_rs.area .+ floe.area)
            for colname in properties
                df_rs[!, Symbol("s1_", colname)] .= floe[colname]
            end
            push!(results, df_rs)
        end
    end
    if length(results) == 0
        return DataFrame(Dict(x=>[] for x in union(properties, [:s1_label, :s2_label, :dist_s1_s2, :scaled_relative_error_area])))
    end
    results_df = vcat(results...; cols=:union)
    rename!(results_df, Dict(r => Symbol("s2_", r) for r in properties))

    return results_df
end

"""
    merge_floes(df1, df2, labels1, labels2)

Produce a single segmentation from a pair via object-wise assessment.
1. Where the two segmentations agree within tolerance of dmax, emax, select the most circular floe.
2. Where the segmentations disagree, select floes with the highest boundary contrast within their

"""
function merge_floes(df1, df2, labels1, labels2; 
    max_distance_pixels=10,
    max_error_area=0.25,
    min_floe_size=100
    )

    # If no floes to merge, skip merge
    nrow(df1) == 0 && return labels2
    nrow(df2) == 0 && return labels1

    #### Set up starting images
    A = labels1
    B = labels2
    offset_b = maximum(A) # Offset the labels in B by the largest value in A
    A_indices = component_indices(A)
    B_indices = component_indices(B)
    A_labels = df1.label
    B_labels = df2.label

    F = zeros(Int64, size(A))

    #### Case 1: No overlap
    A_no_overlap = []
    B_no_overlap = []
    for L in A_labels
        if maximum(B[A_indices[L]]) == 0
            F[A_indices[L]] .= L
            push!(A_no_overlap, L)
        end
    end
    for L in B_labels
        if maximum(A[B_indices[L]]) == 0
            F[B_indices[L]] .= L + offset_b
            push!(B_no_overlap, L)
        end
    end

    subset!(df1, :label => ByRow(r -> r ∉ A_no_overlap))
    subset!(df2, :label => ByRow(r -> r ∉ B_no_overlap))
    nrow(df1) == 0 || nrow(df2) == 0 && return F

    #### Case 2: High-Quality Pairs
    # In this case, there exists at least one item in the relevant set where the error metrics are both within the tolerance.
    # Out of these objects, choose the one with the highest probability. 
    df_comp = objectwise_compare_segmentation(df1, df2, labels1, labels2);
    within_tolerance(d, e) = (d .< max_distance_pixels) .&& (e .< max_error_area)
    matches = subset(df_comp, [:dist_s1_s2, :scaled_relative_error_area] => within_tolerance)
    
    nrow(matches) > 0 && begin
        # Select the item in the relative set with lowest area difference.
        subset!(
            groupby(matches, :s1_label),
            :scaled_relative_error_area => r -> 1:length(r) .== argmin(r),
        )
        subset!(
            groupby(matches, :s2_label),
            :scaled_relative_error_area => r -> 1:length(r) .== argmin(r),
        )

        # Select the option with highest probability
        transform!(
            matches,
            [:s1_probability, :s2_probability] =>
                ByRow((s1, s2) -> s1 .> s2) => :s1_better,
        )

        # Merge the two, prioritizing the second if there is overlap.
        A_labels = matches[matches.s1_better, :s1_label]
        B_labels = matches[.!matches.s1_better, :s2_label];

        for L in A_labels
            F[A_indices[L]] .= L
        end        
        for L in B_labels
            F[B_indices[L]] .= L + offset_b
        end

        # Add intersections to list
        idx = F .> 0
        A_labels = union(A_labels, unique(A[idx]))
        B_labels = union(B_labels, unique(B[idx]))

        # Update the dataframes to remove the resolved labels
        subset!(df1, :label => ByRow(r -> r ∉ A_labels))
        subset!(df2, :label => ByRow(r -> r ∉ B_labels))
    end

    #### Case 3: Poor matches, including over and undersegmentation
    # 1. Loop through remaining objects in A. If probability is higher
    #    for the object in A than all intersections in B, keep object.
    # 2. Loop through remaining objects in B. If no intersection with
    #    the objects kept in step 1, keep object.
    # 3. Update F and return.

    # Select objects in A with higher probability than any intersection with B
    A_labels = []
    B_probability = Dict(r => p for (r, p) in zip(df2.label, df2.probability))
    for s1 in eachrow(df1)        
        B_labels = filter(r -> r ∈ df2.label, unique(labels2[A_indices[s1.label]]))
        if all(s1.probability .> [B_probability[r] for r in B_labels])
            push!(A_labels, s1.label)
        end
    end
    for L in A_labels
        F[A_indices[L]] .= L
    end
    
    # Select objects in B with no intersection with F
    B_labels = unique(B[F .> 0])
    subset!(df2, :label => ByRow(r -> r ∉ B_labels))
    for L in df2.label
        F[B_indices[L]] .= L + offset_b
    end
    
    # Remove possible isolated pixels from merge
    remove_small_segments!(F, min_floe_size)
    return label_components(F)
end

#### Tracker parameters ####
# Initially same as in src/Tracking/filter_functions.jl,
# parameters and filters will be updated based on calibration
# tests.

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
