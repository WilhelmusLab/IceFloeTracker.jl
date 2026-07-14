module FSPipeline
"""

Simplified segmentation pipeline with calibrated parameters for the Greenland Sea / Fram Strait workflow.

"""

using Images
import Dates: Day
import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering: 
    nonlinear_diffusion, 
    PeronaMalikDiffusion, 
    unsharp_mask, 
    ContrastLimitedAdaptiveHistogramEqualization

import ..Morphology:
    fill_holes,
    strel_disk

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
    kmeans_binarization,
    tiled_adaptive_binarization,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionBrightnessMidpoint,
    remove_small_segments!,
    remove_large_segments!,
    stitch_clusters,
    view_seg,
    view_seg_random

import ..Tracking:
    ChainedFilterFunction,
    DistanceThresholdFilter,
    FloeTracker,
    LogLogQuadraticTimeDistanceFunction,
    MinimumWeightMatchingFunction,
    PiecewiseLinearThresholdFunction,
    RelativeErrorThresholdFilter,
    ShapeDifferenceThresholdFilter,
    PsiSCorrelationThresholdFilter

import ..Pipeline: IceFloeSegmentationAlgorithm

abstract type IceFloePreprocessingAlgorithm end

# TODO: In calibration, find a minimum usable rblock/cblock size in pixels.
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
    diffusion_algorithm = PeronaMalikDiffusion(λ=0.1, K=0.1, niters=5, g="exponential")
    adapthisteq_params = (nbins=256, rblocks=8, cblocks=4, clip=1)
    unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
end

function (p::Preprocess)(
    truecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}}, 
    landmask,
    tiles
)
    # Cast to grayscale first to save compute time
    proc_img = Gray.(truecolor_image)
    
    # Diffusion and sharpening
    nonlinear_diffusion(proc_img, tiles, p.diffusion_algorithm)
   
    adjust_histogram!(proc_img,
        ContrastLimitedAdaptiveHistogramEqualization(
            nbins=p.adapthisteq_params.nbins,
            rblocks=p.adapthisteq_params.rblocks,
            cblocks=p.adapthisteq_params.cblocks,
            clip=p.adapthisteq_params.clip)
    )

    proc_img .= unsharp_mask(proc_img,
        p.unsharp_mask_params.radius,
        p.unsharp_mask_params.amount,
        p.unsharp_mask_params.threshold
    )
            
    # Re-apply mask so sharpening doesn't bleed into land
    apply_landmask!(proc_img, landmask)
    return proc_img
end

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
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))
    cloud_mask_algorithm = Watkins2026CloudMask()
    preprocessing_algorithm = Preprocess(
        diffusion_algorithm = PeronaMalikDiffusion(λ=0.1, K=0.1, niters=5, g="exponential"),
        adapthisteq_params = (nbins=256, rblocks=8, cblocks=4, clip=1),
        unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01),
    )
    tile_size_pixels = 1200
    min_tile_ice_pixel_count=300
    preliminary_ice_mask = IceDetectionBrightnessMidpoint(
        minimum_reflectance=0.3
        )
    kmeans_params = (
        k=4, 
        maxiter=50, 
        random_seed=45, 
        cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
            band_7_max=0.1,
            possible_ice_threshold=0.3,
            join_method="union",
            minimum_prominence=0.01,
            )
        )
    adaptive_params = (
        window_size=400,
        percentage=0
    )
    cleanup_binary_params = (
        erosion_strel=strel_box((3,3)), 
        init_max_fill=100, 
        conditional_max_fill=500
    )
    floe_splitting_params = (
        max_hole_fill=2000,
        max_distance=5,
        max_expand=3,
    )
    floe_filtering_params = (
        min_floe_size=100,
        max_floe_size=50_000,
        min_band_2_reflectance=0.4,
        max_band_7_reflectance=0.7,
        min_circularity=0.4
    )
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
    apply_landmask!(truecolor_image, landmask)
    apply_landmask!(falsecolor_image, landmask)

    n, m = size(truecolor_image)
    tile_size_pixels = s.tile_size_pixels
    tile_size_pixels > maximum([n, m]) && begin
        @warn "Tile size too large, defaulting to image size"
            tile_size_pixels = minimum([n, m])
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
                t -> sum(.!joint_mask[t...]) > s.min_tile_ice_pixel_count, tiles);

    # Then check for sufficient possible sea ice pixels
    prelim_ice_mask = s.preliminary_ice_mask(Gray.(red.(tc_masked)), filtered_tiles)
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > s.min_tile_ice_pixel_count, filtered_tiles);

    @info "Preprocessing truecolor image"
    preproc_gray = float64.(s.preprocessing_algorithm(
        truecolor_image, landmask, filtered_tiles));

    @info "Binarization"
    # We use the cloud mask in finding the bright floes - the bright floe cluster can't be cloud -
    # and allow the k-means cluster to overlap with the cloud mask by using the preproc gray with
    # only the landmask applied to it
    # Update to do the k-means twice: preproc and not preproc
    kmeans_result = kmeans_binarization(
            apply_landmask(preproc_gray, cloud_mask),
            fc_masked,
            filtered_tiles;
            s.kmeans_params...
            )
    adaptive_result = binarize(preproc_gray,
        AdaptiveThreshold(;
        s.adaptive_params... 
        )
    ) .> 0
    
    # AdaptiveThreshold often has noise in large blank areas
    apply_landmask!(adaptive_result, landmask)


    @info "Splitting floes"

    clean_split_label = r -> dist_morph_split(
                clean_binary_floes(r, prelim_ice_mask, cloud_mask; s.cleanup_binary_params...
            );
            s.floe_splitting_params...
        )
    
    kmeans_split_floes = clean_split_label(kmeans_result)
    adaptive_split_floes = clean_split_label(adaptive_result)
    
    # TBD: Filter floes based on the edge properties, colors

    @info "Filtering floes"
    
    # Includes removal of objects intersecting the coastal buffer and
    # objects outside the range (min_floe_size, max_floe_size)
    # TBD: object-based analysis

    filter_floes_! = r -> filter_floes!(
        r,
        coastal_buffer_mask, 
        falsecolor_image,
        s.floe_filtering_params.min_floe_size, 
        s.floe_filtering_params.max_floe_size,
        s.floe_filtering_params.min_band_2_reflectance,
        s.floe_filtering_params.max_band_7_reflectance,
        s.floe_filtering_params.min_circularity,
    )

    filter_floes_!(kmeans_split_floes)
    filter_floes_!(adaptive_split_floes)
   

    # Object-wise analysis
    # - Circularity measure
    # - Border analysis

    # Join results
    final_floes = kmeans_split_floes

    # Re-label so there are no missing numbers in the component list
    final_floes .= label_components(final_floes)

    # Return the original truecolor image, segmented
    segments_tc = SegmentedImage(truecolor_image, final_floes)
    segments_fc = SegmentedImage(falsecolor_image, final_floes)

    if !isnothing(intermediate_results_callback)
        colorview_truecolor = view_seg(segments_tc)
        colorview_falsecolor = view_seg(segments_fc)
        colorview_random =  view_seg_random(segments_tc)
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            coastal_buffer_mask=Gray.(coastal_buffer_mask),
            cloud_mask=Gray.(cloud_mask),
            ice_mask=Gray.(prelim_ice_mask),
            preprocessed=preproc_gray,
            kmeans_binarized=kmeans_split_floes .> 0,
            adaptive_binarized=adaptive_split_floes .> 0,
            final_floes = colorview_random,
            labels_map = final_floes,
            segment_mean_falsecolor=colorview_falsecolor,
            segment_mean_truecolor=colorview_truecolor,
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
function clean_binary_floes(binary_img, icemask, cloudmask;
        erosion_strel=strel_box((3,3)),
        init_max_fill=100,
        conditional_max_fill=500
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
        opening_strel=strel_disk(3)
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
            (L > 0) && begin
                matched_labels = unique(levels[dist_threshold][indices[L]])
                
                # If intersection of the label at level 
                if (0 ∈ matched_labels) && (length(matched_labels) <= 2)
                    final_labels[indices[L]] .= L
                else
                    # Otherwise, expand the current level, and set the next level down to the expanded indices.
                    # May need to check the number of matched labels in the expanded image.
                    levels[dist_threshold - 1][indices[L]] .= expanded[indices[L]]
                    final_labels[indices[L]] .= expanded[indices[L]]
                end
            end
        end
    end
    return label_components(final_labels)
end

# Helper function for creating a filtered version of the image indexmap
"""assign_labels(img_indexmap, labels_list)

Given an image indexmap `img` and a `labels_list`, create a new labeled
image using only the values in the list.

"""
function assign_labels(img_indexmap, labels_list)
    out = zeros(Int64, size(img_indexmap))
    indices = component_indices(img_indexmap)
    for L in labels_list
        out[indices[L]] .= L
    end
    return out
end  


function filter_floes!(
    img_indexmap, 
    coastal_buffer_mask, 
    falsecolor_image,
    min_floe_size, 
    max_floe_size,
    min_band_2_reflectance,
    max_band_7_reflectance,
    min_circularity,
    )
    overlap = unique(img_indexmap[coastal_buffer_mask])
    indices = component_indices(img_indexmap)    
    for L in overlap
        img_indexmap[indices[L]] .= 0        
    end

    # Remove floes outside the specified bounds
    remove_small_segments!(img_indexmap, min_floe_size)
    remove_large_segments!(img_indexmap, max_floe_size)

    areas = component_lengths(img_indexmap)
    perims = component_perimeters(img_indexmap)
    b2_means = segment_mean(SegmentedImage(green.(falsecolor_image), img_indexmap))
    b7_means = segment_mean(SegmentedImage(red.(falsecolor_image), img_indexmap))

    for L in unique(img_indexmap)
        if L > 0 
            if b2_means[L] < min_band_2_reflectance || b7_means[L] > max_band_7_reflectance
                img_indexmap[indices[L]] .= 0
            elseif 4 * π * areas[L] / perims[L]^2 < min_circularity
                img_indexmap[indices[L]] .= 0
            end
        end
    end
end

# TODO: Filter function for object-wise cleanup
# TODO: Add object-wise hole filling method
# TODO: Update settings for tracking based on test cases

# TODO: Add this to the segmentation utilities
"""
    get_relevant_set(ground_truth_df, segments_df, ground_truth_labels, segment_labels)   


Find the relevant set for comparing ground truth labels with a segmentation product.
- ground_truth_df = result of regionprops table
- segments_df = result of regionprops table
- ground_truth_labels = labels map for the ground truth
- segment_labels = labels map for the segmentation results

The relevant set for a segmentation comparison is defined by 
1. centroid g in s
2. centroid s in g
3. area overlap greater than 50% of g
4. area overlap greater than 50% of s

This method focuses on floes which overlap, not false positives.

"""
function get_relevant_set(ground_truth_df, segments_df, ground_truth_labels, segment_labels)   
    relevant_set = Dict{Int64, Vector{Int64}}()  
    for floe in eachrow(ground_truth_df)
        # select labels that are inside the bounding box for the floe
        matched_labels = unique(segment_labels[floe.min_row:floe.max_row, floe.min_col:floe.max_col])

        # if any, then check centroid positions
        maximum(matched_labels) != 0 && begin
            # get the rows in the segments_df from the matched labels
            candidate_subset = subset(segments_df, :label => ByRow(l -> l in matched_labels))

            relevant_set_labels = []
            
            # check if centroid g in s
            rc = round(Int64, floe.row_centroid)
            cc = round(Int64, floe.col_centroid)
            push!(relevant_set_labels, segment_labels[rc, cc])

            # check if centroid s in g
            for s_floe in eachrow(candidate_subset)
                rc = round(Int64, s_floe.row_centroid)
                cc = round(Int64, s_floe.col_centroid)
                (ground_truth_labels[rc, cc] == floe.label) && begin
                    push!(relevant_set_labels, s_floe.label)
                end

                # joint bbox
                rmin = minimum((floe.min_row, s_floe.min_row))
                rmax = maximum((floe.max_row, s_floe.max_row))
                cmin = minimum((floe.min_col, s_floe.min_col))
                cmax = maximum((floe.max_col, s_floe.max_col))
                
                # check if area overlap between g and s is larger than 50% of g
                gtmask = ground_truth_labels[rmin:rmax, cmin:cmax] .== floe.label
                slmask = segment_labels[rmin:rmax, cmin:cmax] .== s_floe.label
                intersect_area = sum(gtmask .&& slmask)
                maximum([intersect_area / s_floe.area, intersect_area / floe.area]) > 0.5 && begin
                    push!(relevant_set_labels, s_floe.label)
                end
            end
            relevant_set_labels = filter(r -> r != 0, unique(relevant_set_labels))
            if length(relevant_set_labels) > 0
                push!(relevant_set, floe.label => relevant_set_labels)
            end
        end
    end
    return relevant_set
end

"""
Objectwise comparison: Get the index for the relevant set,
compute the centroid distance, relative error in area, and 
add info on the matched segment.
"""
function objectwise_compare_segmentation(segmentation1, segmentation2)

    # Get the region properties 
    # Extend this as we add more measures
    # - Shape difference
    # - Intersection area
    # - Boundary distance
    
    df_s1 = regionprops_table(segmentation1; properties=[:label, :centroid, :area, :bbox, :perimeter])
    df_s2 = regionprops_table(segmentation2; properties=[:label, :centroid, :area, :bbox, :perimeter])

    relevant_set = get_relevant_set(df_s1, df_s2, labels_map(segmentation1), labels_map(segmentation2))
    results = DataFrame[]
    for floe in eachrow(df_s1)
        g = floe.label
        g in keys(relevant_set) && begin
            df_rs = subset(df_s2, :label => ByRow(s -> s in relevant_set[g]))        
            df_rs[:, :s1_label] .= g
            df_rs[:, :s1_area] .= floe.area
            df_rs[:, :s1_perimeter] .= floe.perimeter
            df_rs[:, :s1_row_centroid] .= floe.row_centroid
            df_rs[:, :s1_col_centroid] .= floe.col_centroid
            df_rs[:, :dist_s1_s2] = euclidean_distance(floe, df_rs)
            df_rs[:, :scaled_relative_error_area] = abs.(df_rs.area .- floe.area) ./ (df_rs.area .+ floe.area)
            push!(results, df_rs)
        end
    end
    results_df = vcat(results..., cols=:union)
    rename!(results_df, :area => :s2_area, :perimeter => :s2_perimeter,  :label => :s2_label, :col_centroid => :s2_col_centroid, 
                        :row_centroid => :s2_row_centroid, :max_col => :s2_max_col,
                        :max_row => :s2_max_row, :min_col=>:s2_min_col, :min_row=>:s2_min_row)

    @. results_df[:, :s1_circularity] = 4 * pi * results_df[:, :s1_area] / results_df[:, :s1_perimeter]^ 2
    @. results_df[:, :s2_circularity] = 4 * pi * results_df[:, :s2_area] / results_df[:, :s2_perimeter]^ 2
    
    return_cols = [
         "s1_label",
         "s1_area",
         "s1_perimeter",
         "s1_row_centroid",
         "s1_col_centroid",
         "s1_circularity",
         "s2_label",
         "s2_area",
         "s2_perimeter", 
         "s2_col_centroid",
         "s2_row_centroid",
         "s2_circularity",
         "dist_s1_s2",
         "scaled_relative_error_area"
    ]

    return results_df[:, return_cols]
end



"""
    Track()

Track shapes across images using the LogLogQuadratic distance filter, the ChainedFilterFunction,
and the MinimumWeightMatchingFunction.

"""
function Track(
    filter_function=ChainedFilterFunction(;
        filters=[
            DistanceThresholdFilter(
                threshold_function=LogLogQuadraticTimeDistanceFunction(),
            ),
            RelativeErrorThresholdFilter(;
                variable=:area,
                threshold_function=PiecewiseLinearThresholdFunction(;
                    minimum_area = 100,
                    maximum_area = 700,
                    minimum_value=0.43,
                    maximum_value=0.17,
                ),
            ),
            RelativeErrorThresholdFilter(;
                variable=:convex_area,
                threshold_function=PiecewiseLinearThresholdFunction(;
                    minimum_area = 100,
                    maximum_area = 700,
                    minimum_value=0.44,
                    maximum_value=0.25,
                ),
            ),
            RelativeErrorThresholdFilter(;
                variable=:major_axis_length,
                threshold_function=PiecewiseLinearThresholdFunction(;
                    minimum_area = 100,
                    maximum_area = 700,
                    minimum_value=0.27,
                    maximum_value=0.13,
                ),
            ),
            RelativeErrorThresholdFilter(;
                variable=:minor_axis_length,
                threshold_function=PiecewiseLinearThresholdFunction(;
                    minimum_area = 100,
                    maximum_area = 700,
                    minimum_value=0.28,
                    maximum_value=0.1,
                ),
            ),
            ShapeDifferenceThresholdFilter(;
                threshold_function=PiecewiseLinearThresholdFunction(;
                    minimum_area = 100,
                    maximum_area = 700,
                    minimum_value=0.47,
                    maximum_value=0.31,
                ),
            ),
            PsiSCorrelationThresholdFilter(;
                threshold_function=PiecewiseLinearThresholdFunction(;
                    minimum_area = 100,
                    maximum_area = 700,
                    minimum_value=0.86,
                    maximum_value=0.96,
                ),
            ),
        ],
    ),
    matching_function=MinimumWeightMatchingFunction(
        columns = [
        :scaled_distance,
        :relative_error_area,
        :relative_error_convex_area,
        :relative_error_major_axis_length,
        :relative_error_minor_axis_length,
        :psi_s_correlation_score,
        :scaled_shape_difference,
        ],
        weights = ones(7)
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
