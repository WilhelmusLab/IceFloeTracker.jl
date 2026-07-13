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
    LopezAcostaCloudMask,
    Watkins2025CloudMask
import ..ImageUtils: get_tiles, imbrighten
import ..Segmentation: 
    expand_labels,
    kmeans_binarization,
    tiled_adaptive_binarization,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionBrightnessMidpoint,
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
    adapthisteq_params = (nbins=256, rblocks=8, cblocks=4, clip=5)
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
- `min_floe_size=100`: Smallest floe size to retain
- `preliminary_ice_mask = IceDetectionBrightnessMidpoint(minimum_reflectance=0.3)`: Function to use to identify likely ice pixels for filtering.
- `kmeans_params = (k=4, maxiter=50, random_seed=45)`: Parameters for `kmeans_binarization`
- `cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
    band_7_max=0.1,
    possible_ice_threshold=0.3,
    join_method="union",
    minimum_prominence=0.01)`: Function to use to select a k-means cluster in the `kmeans_binarization` workflow
- `clean_binary_floes_params`: Parameters for the preliminary binary image cleanup
- `floe_splitting_params`: Parameters for the `dist_morph_split` floe splitting algorithm
"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm 
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))
    cloud_mask_algorithm = Watkins2025CloudMask()
    preprocessing_algorithm = Preprocess()
    tile_size_pixels = 1200
    min_tile_ice_pixel_count=300
    min_floe_size=100
    max_floe_size=50_000
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    preliminary_ice_mask = IceDetectionBrightnessMidpoint(minimum_reflectance=0.3)
    cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01)
    floe_splitting_settings = (max_fill_area=1, min_area_opening=20, opening_strel=strel_disk(2))
    # TBD: Add updated floe splitting settings, add params for k-means cleanup
end 


function (p::Segment)(
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
    tile_size_pixels = p.tile_size_pixels
    tile_size_pixels > maximum([n, m]) && begin
        @warn "Tile size too large, defaulting to image size"
            tile_size_pixels = minimum([n, m])
    end
    
    (nr, nc) = round.(Int, size(truecolor_image) ./ tile_size_pixels)
    tiles = get_tiles(truecolor_image; rblocks=nr, cblocks=nc)
    @info "Building masks"
    cloud_mask = create_cloudmask(falsecolor_image, p.cloud_mask_algorithm)
    
    # 2. Intermediate images - apply coastal buffer and cloud mask
    joint_mask = coastal_buffer_mask .|| cloud_mask
    tc_masked = apply_landmask(truecolor_image, joint_mask)
    fc_masked = apply_landmask(falsecolor_image, joint_mask)

    # First check for sufficient non-land and non-cloud pixels
    filtered_tiles = filter(
                t -> sum(.!joint_mask[t...]) > p.min_tile_ice_pixel_count, tiles);

    # Then check for sufficient possible sea ice pixels
    prelim_ice_mask = p.preliminary_ice_mask(tc_masked, filtered_tiles)
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > p.min_tile_ice_pixel_count, filtered_tiles);

    @info "Preprocessing truecolor image"
    preproc_gray = float64.(p.preprocessing_algorithm(
        truecolor_image, landmask, filtered_tiles));

    # We use the cloud mask in finding the bright floes - the bright floe cluster can't be cloud -
    # and allow the k-means cluster to overlap with the cloud mask by using the preproc gray with
    # only the landmask applied to it
    # Alternative approach: simply use the adaptive threshold binarization. 
    kmeans_result = kmeans_binarization(
            apply_landmask(preproc_gray, cloud_mask),
            fc_masked,
            filtered_tiles;
            k=p.kmeans_params.k,
            maxiter=p.kmeans_params.maxiter,
            random_seed=p.kmeans_params.random_seed,
            cluster_selection_algorithm=p.cluster_selection_algorithm
            )
     # update to have settings accessible from top
    kmeans_result .= clean_binary_floes(kmeans_result, prelim_ice_mask, cloud_mask)
    # kmeans_result = tiled_adaptive_binarization(apply_landmask(preproc_gray, cloud_mask),
    #     filtered_tiles;
    #     minimum_window_size=400,
    #     minimum_brightness=0.3,
    #     threshold_percentage=0
    #     ) .> 0 
    kmeans_result .= clean_binary_floes(kmeans_result, prelim_ice_mask, cloud_mask)

    @info "Splitting floes"
    # Could tile this, but doesn't seem to be a major bottleneck
    split_floes = dist_morph_split(kmeans_result; max_distance=7) # update to have morph split settings
    # TBD: Filter floes based on the edge properties, colors

    @info "Filtering floes"
    
    # Remove floes which intersect the coastal buffer
    # This one could be a separate function like remove_small_segments!
    overlap = unique(split_floes[coastal_buffer_mask])
    indices = component_indices(split_floes)    
    for L in overlap
        split_floes[indices[L]] .= 0        
    end

    remove_small_segments!(split_floes, p.min_floe_size)
    remove_large_segments!(split_floes, p.max_floe_size)
    # Re-label to fill where regions were deleted
    split_floes .= label_components(split_floes)

    # Return the original truecolor image, segmented
    segments_tc = SegmentedImage(truecolor_image, split_floes)
    segments_fc = SegmentedImage(falsecolor_image, split_floes)

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
            bright_ice_mask=p.cluster_selection_algorithm(falsecolor_image),            
            binarized=kmeans_result,
            final_floes = colorview_random,
            labels_map = split_floes,
            segment_mean_falsecolor=colorview_falsecolor,
            segment_mean_truecolor=colorview_truecolor,
            ) 
    end
    return segments_tc
end

function clean_binary_floes(binary_img, icemask, cloudmask;
        erosion_strel=strel_box((7,7)),
        filling_strel=strel_diamond((3,3)),
        max_fill=100
    )
    out = deepcopy(binary_img)
    # 1. Shrink objects using the provided structuring element
    eroded_img = erode(out, erosion_strel)

    # 2. After shrinking, fill holes
    filled = fill_holes(eroded_img, filling_strel) # Test how permissive this is. Should we use imfill instead?

    # 3. Identify filled holes which are part of the ice mask or the cloud mask
    filled .= filled .&& (icemask .|| cloudmask)
    filled .= .!imfill(.!filled, (0, max_fill))

    # 4. Use morphological closing to further limit openings
    closing!(filled)

    # 5. Finally, set any of these filled pixels to 1 in the output image.
    out[filled .> 0] .= 1
    return out
end


### TODO: Set up this function as a "FloeSplittingAlgorithm" 
# Find markers by selecting locations greater than dist threshold from background
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
        min_floe_size::Int64=64, # TBD: add maximum floe size
        max_hole_fill::Int64=2000,
        max_distance::Int64=5,
        max_expand::Int64=3,
        opening_strel=strel_disk(3)
    )

    ### Remove objects below the minimum floe size, then create distance pyramid
    bw = .!imfill(binary_floes, (0, min_floe_size))
    dist = distance_transform(feature_transform(bw))
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
    final_labels .= label_components(final_labels)
    remove_small_segments!(final_labels, min_floe_size)
    return final_labels
end

# TODO: Add function to segmented image utilities and add test
"""
    remove_small_segments!(labels, min_size)

Checks the area of each labeled object in `labels` and sets it to 0 if it is less than `min_size`.

"""
function remove_small_segments!(labels, min_size)
    areas = component_lengths(labels)
    indices = component_indices(labels)

    for L in keys(areas)
        (L != 0) && begin
            (areas[L] < min_size) && begin
                labels[indices[L]] .= 0
            end
        end
    end
end

"""
    remove_small_segments!(labels, min_size)

Checks the area of each labeled object in `labels` and sets it to 0 if it is less than `min_size`.

"""
function remove_large_segments!(labels, max_size)
    areas = component_lengths(labels)
    indices = component_indices(labels)

    for L in keys(areas)
        (L != 0) && begin
            (areas[L] > max_size) && begin
                labels[indices[L]] .= 0
            end
        end
    end
end

#### TODO: Add object-wise hole filling method
#### TODO: Add Track() method for configured tracker (including alternative similarity measures)

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
