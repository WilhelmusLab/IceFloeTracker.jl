module FSPipeline
"""
Update of the LopezAcosta2019 algorithm for the Fram Strait floe trajectory dataset. Key updates:
    * Separating segmentation algorithm components for clearer functions
    * Implementation of the tiling workflow
    * Simpler watershed transformation steps
    * Exposure of algorithm settings to the user
"""

using Images
import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering: 
    nonlinear_diffusion, 
    PeronaMalikDiffusion, 
    unsharp_mask, 
    ContrastLimitedAdaptiveHistogramEqualization

import ..Morphology: hbreak, hbreak!, branch, bridge, fill_holes, strel_disk
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
    IceFloeSegmentationAlgorithm, 
    find_ice_mask, 
    kmeans_binarization,
    tiled_adaptive_binarization,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionBrightnessPeaksMODIS134,
    stitch_clusters,
    view_seg

abstract type IceFloePreprocessingAlgorithm end

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

"""
@kwdef struct Preprocess <: IceFloePreprocessingAlgorithm
    diffusion_algorithm = PeronaMalikDiffusion(λ=0.1, K=0.1, niters=5, g="exponential")
    adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.99) # rblocks/cblocks not used yet -- add with CLAHE.jl
    unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
end

function (p::Preprocess)(
    truecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}}, 
    landmask,
    cloud_mask
)
    # Cast to grayscale first to save compute time
    proc_img = Gray.(truecolor_image)
    apply_landmask!(proc_img, landmask .|| cloud_mask)

    # Diffusion and sharpening
    proc_img .= nonlinear_diffusion(proc_img, p.diffusion_algorithm)

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
            
    apply_landmask!(proc_img, landmask .|| cloud_mask)
    return proc_img
end


"""
    FSPipeline.Segment()

Segmentation routine for identifying moderate to large floes in the Fram Strait. 
The image preprocessing is supplied as an function in the functor setup.


# Parameters
- `coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))`: Structuring element for the `create_landmask` function
- `cloud_mask_algorithm = LopezAcostaCloudMask()`: Cloud mask algorithm
- `preprocessing_algorithm = Preprocess()`: Function to sharpen and equalize the truecolor image
- `tile_size_pixels=1000`: Nominal tile size in pixels
- `min_tile_ice_pixel_count=300`: Smallest number of likely sea ice floes in tile
- `min_floe_size=100`: Smallest floe to retain in image
- `preliminary_ice_mask = IceDetectionBrightnessPeaksMODIS134(band_7_max=0.1, possible_ice_threshold=0.3)`: Function to use to identify likely ice pixels for filtering.
- `kmeans_params = (k=4, maxiter=50, random_seed=45)`: Parameters for `kmeans_binarization`
- `cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
    band_7_max=0.1,
    possible_ice_threshold=0.3,
    join_method="union",
    minimum_prominence=0.01)`: Function to use to select a k-means cluster in the `kmeans_binarization` workflow
- `brightening_factor = 0.3`
- `adaptive_binarization_settings = (
    minimum_window_size=400,
    threshold_percentage=0,
    minimum_brightness=100/255)`: Settings for the adaptive threshold binarization
- `watershed_strel = se_disk(5)`
- `floe_splitting_settings = (max_fill_area=1, min_area_opening=20, opening_strel=se_disk(2))`
"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm # Tried making this an IceFloeSegmentationAlgorithm but julia complains about ambiguity when I do so. Need to fix this to be able to use the run and validate function.
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))
    cloud_mask_algorithm = Watkins2025CloudMask()
    preprocessing_algorithm = Preprocess()
    tile_size_pixels=1200
    min_tile_ice_pixel_count=300
    min_floe_size=100
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    preliminary_ice_mask = IceDetectionBrightnessPeaksMODIS134(band_1_min=0.3)
    cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01)
    floe_splitting_settings = (max_fill_area=1, min_area_opening=20, opening_strel=strel_disk(2))
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

    n, m = size(truecolor_image)
    tile_size_pixels = p.tile_size_pixels
    tile_size_pixels > maximum([n, m]) && begin
        @warn "Tile size too large, defaulting to image size"
            tile_size_pixels = minimum([n, m])
    end
    # TODO: Option to supply number of row and cols 
    tiles = get_tiles(truecolor_image, tile_size_pixels)
 
    @info "Building masks"
    # TODO: Make sure tests aren't over-sensitive to roundoff errors for Float32 vs Float64
    cloud_mask = create_cloudmask(falsecolor_image, p.cloud_mask_algorithm)
    landmask = landmask .> 0 # make sure it's a bitmatrix
    # 2. Intermediate images - using coastal buffer on the FC image
    apply_landmask!(truecolor_image, landmask)
    apply_landmask!(falsecolor_image, coastal_buffer_mask .|| cloud_mask)

    prelim_ice_mask = p.preliminary_ice_mask(truecolor_image, tiles)
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > p.min_tile_ice_pixel_count, tiles);

    @info "Preprocessing truecolor image"
    preproc_gray = float64.(p.preprocessing_algorithm(
        truecolor_image, landmask, landmask, filtered_tiles));
        # Uses the landmask twice here -- intentionaly bypassing the cloud mask!
    
    # We use the cloud mask in finding the bright floes - the bright floe cluster can't be cloud.
    kmeans_result = kmeans_binarization(
            preproc_gray,
            falsecolor_image, 
            filtered_tiles;
            k=p.kmeans_params.k,
            maxiter=p.kmeans_params.maxiter,
            random_seed=p.kmeans_params.random_seed,
            cluster_selection_algorithm=p.cluster_selection_algorithm
            )

    kmeans_result .= clean_binary_floes(kmeans_result, prelim_ice_mask, cloud_mask) # update to have settings

    @info "Splitting floes"
    # how to tile this?
    split_floes = dist_morph_split(kmeans_result; max_distance=7) # update to have morph split settings
    
    # Remove too-small floes

    # Filter floes based on the edge properties, colors

    # re-label 
    labels = label_components(labels) # reset labels

    # Return the original truecolor image, segmented
    segments = SegmentedImage(truecolor, labels)

    if !isnothing(intermediate_results_callback)
        colorview_truecolor = view_seg(SegmentedImage(truecolor, labels))
        colorview_falsecolor = view_seg(SegmentedImage(falsecolor, labels))
        colorview_random =  view_seg_random(SegmentedImage(truecolor, labels))
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
            segment_mean_falsecolor=colorview_falsecolor,
            segment_mean_truecolor=colorview_truecolor,
            ) 
    end
    return segments
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

# Expand labels by distance without overlap
function expand_labels(labels, distance)
    labels_out = deepcopy(labels)
    maximum(labels_out) == 0 && return labels_out
    F = feature_transform(labels .> 0)
    D = distance_transform(F)
    labels_out[D .<= distance] .= labels[F][D .<= distance]
    return labels_out
end

# Find markers by selecting locations greater than dist threshold from background
function dist_morph_split(
        binary_floes::BitMatrix;
        min_floe_size::Int64=64,
        max_hole_fill::Int64=2000,
        max_distance::Int64=5,
        max_expand::Int64=3,
        strel=se_disk(3)
    )
    bw = .!imfill(binary_floes, (0, min_floe_size))
    dist = distance_transform(feature_transform(bw))
    levels = Dict(0 => label_components(opening(dist .> 0, strel))) # Initialize with one run of opening
    ### Build pyramid ###
    for dist_threshold in 1:max_distance
        markers = opening(dist .> dist_threshold, strel)
        markers .= .!imfill(.!markers, (0, 2000))
        levels[dist_threshold] = label_components(markers)
    end
    final_labels = deepcopy(levels[max_distance])

    ### Descend pyramid ####
    for dist_threshold in max_distance:-1:1
        # Get indices of next level down
        indices = component_indices(levels[dist_threshold - 1])

        # Expand indices of current level
        expanded = expand_labels(levels[dist_threshold], max_expand)
        for L in keys(indices)
            (L > 0) && begin
                matched_labels = unique(levels[dist_threshold][indices[L]])
                
                # If no higher levels or only one higher level, set to current label
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
    areas = component_lengths(final_labels)
    indices = component_indices(final_labels)
    for L in keys(areas)
        (areas[L] < min_floe_size) && (final_labels[indices[L]] .= 0)
    end
    return final_labels
end

end
