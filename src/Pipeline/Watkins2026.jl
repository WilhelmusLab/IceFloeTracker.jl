module Watkins2026
using Images # Find individual functions to import later

import ..Filtering: nonlinear_diffusion, PeronaMalikDiffusion, unsharp_mask, channelwise_adapthisteq
import ..Morphology: hbreak, hbreak!, branch, bridge
import ..Preprocessing:
    create_landmask,
    create_cloudmask,
    apply_landmask,
    apply_landmask!,
    apply_cloudmask,
    apply_cloudmask!,
    Watkins2025CloudMask

import ..Segmentation: 
    IceFloeSegmentationAlgorithm, 
    kmeans_binarization,
    tiled_adaptive_binarization,
    IceDetectionBrightnessPeaksMODIS721,
    stitch_clusters,
    get_ice_peaks,
    regionprops,
    component_floes,
    component_convex_areas,
    component_perimeters


import ..ImageUtils:
    get_tiles


abstract type IceFloePreprocessingAlgorithm end

###### Default Parameters ######

###### Preprocessing #####
"""
   Preprocess(
        coastal_buffer_structuring_element = strel_box((51,51))
        cloud_mask_algorithm = Watkins2025CloudMask()
        diffusion_algorithm = PeronaMalikDiffusion(lambda=0.1, kappa=0.1, niters=5, g="exponential")
        adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.95)
        unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
    )
    Preprocess()(img, cloudmask, landmask)

    Converts input image to grayscale, then preprocesses by appling nonlinear diffusion, 
    adaptive histogram equalization, and unsharp masking. Diffusion and unsharp masking are applied 
    to each tile, while the adaptive histogram equalization is divided according to the parameter specifications.

"""
abstract type IceFloePreprocessingAlgorithm end

@kwdef struct Preprocess <: IceFloePreprocessingAlgorithm
    diffusion_algorithm = PeronaMalikDiffusion(λ=0.1, K=0.1, niters=5, g="exponential")
    adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.99) # rblocks/cblocks not used yet -- add with CLAHE.jl
    unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
end

# TODO: Decide how to pass landmasks and cloudmasks around based
function (p::Preprocess)(
    truecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}}, 
    land_mask,
    cloud_mask,
    tiles
)
    proc_img = deepcopy(Gray.(truecolor_image))

    # Diffusion and sharpening
    begin
        for tile in tiles
            proc_img[tile...] .= unsharp_mask(proc_img[tile...],
                    p.unsharp_mask_params.radius,
                    p.unsharp_mask_params.amount,
                    p.unsharp_mask_params.threshold)
            proc_img[tile...] .= nonlinear_diffusion(proc_img[tile...], p.diffusion_algorithm)
        end
    end
    
    apply_cloudmask!(proc_img, cloud_mask)
    apply_landmask!(proc_img, land_mask)

    # Histogram equalization
    # Replace with CLAHE.jl when available
    proc_img = adjust_histogram(proc_img, AdaptiveEqualization(;p.adapthisteq_params...))

    # Reapply masks, since sharpening and eq adjustment can bleed into neighboring pixels
    apply_cloudmask!(proc_img, cloud_mask)
    apply_landmask!(proc_img, land_mask)
    
    return proc_img
end

###### Segmentation ########

@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    preprocessing_function = Preprocess()
    coastal_buffer_structuring_element = strel_box((51,51))
    cloud_mask_algorithm = Watkins2025CloudMask()
    tile_size_pixels = 1000
    min_tile_ice_pixel_count = 4000
    bright_ice_detector = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01
    )
    grayscale_floe_threshold = 0.4 # Brighten pixels brighter than this threshold 
    adaptive_binarization_settings = (minimum_window_size=400, threshold_percentage=0, minimum_brightness=100/255)
    # cleanup parameters: can be named tuple, after writing self-contained cleanup function
    minimum_floe_size = 100
    maximum_fill_size = 64
    #kmeans_binarization_settings = (k=4, )
end

function (s::Segment)(truecolor_image, falsecolor_image, land_mask_image)

    tiles = get_tiles(truecolor_image, s.tile_size_pixels)

    cloud_mask = s.cloud_mask_algorithm(falsecolor_image)
    
    _lm_temp = create_landmask(land_mask_image, s.coastal_buffer_structuring_element)
    land_mask = _lm_temp.non_dilated
    coastal_buffer = _lm_temp.dilated
    
    # Ice-water mask from the red channel reflectance
    # TODO: Add tiles option, and make it an IceDetectionAlgorithm to use as an import
    prelim_ice_mask = ice_water_mask(truecolor_image, cloud_mask, land_mask)

    # Enhanced grayscale image for segmentation
    preprocessed_image = s.preprocessing_function(truecolor_image, land_mask, cloud_mask, tiles)

    # Select a subset of the tiles to preprocess by checking for clear-sky ice pixels outside the 
    # coast buffer
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...] .&& .! coastal_buffer[t...]) > s.min_tile_ice_pixel_count,
        tiles);

    @info "Segmentation method 1: First K-means workflow"
    # TODO: Test impact of the three different "floe brightening" methods
    # - Morphological residue method from LopezAcosta2019Tiling
    # - Discriminate-ice-water from LopezAcosta2019
    # - Discriminate-ice-water from the LASW2019 paper
    # - Gamma threshold method from LopezAcosta2019
    # - Make the kmeans_binarized a binarization algorithm
    begin
        fc_masked = apply_cloudmask(falsecolor_image, cloud_mask)
        # Use the coastal buffer here to highlight bright pixels in the ocean rather than landfast
        apply_landmask!(fc_masked, coastal_buffer)
        kmeans_binarized = kmeans_binarization(preprocessed_image, fc_masked, filtered_tiles;
            cluster_selection_algorithm=s.bright_ice_detector, k=4) # 
        kmeans_binarized .= clean_binary_floes(kmeans_binarized, filtered_tiles)
        kmeans_watershed = watershed_transform(kmeans_binarized, truecolor_image, filtered_tiles; 
                                strel=strel_diamond((5,5)),
                                dist_threshold=4, min_overlap=10,
                                grayscale_threshold=0.3) # output type? I think it's a segmented image                
        kmeans_bdry = falses(size(preprocessed_image))
        for tile in filtered_tiles # Could have a "tiled boundary" function
            kmeans_bdry[tile...] = isboundary((labels_map(kmeans_watershed) .== 0)[tile...]) # Question: boundary of background, or boundary of objects?
        end
    end

    @info "Segmentation method 2: Adaptive threshold workflow" # TODO:make the post processing a function, since it's repeated
    begin
        adaptive_thresh =
                tiled_adaptive_binarization(
                    preprocessed_image,
                    filtered_tiles; s.adaptive_binarization_settings...
                ) .> 0 # returns Gray 
        adaptive_thresh .= clean_binary_floes(adaptive_thresh, filtered_tiles)
        adapt_watershed = watershed_transform(adaptive_thresh, truecolor_image, filtered_tiles;
                                 strel=strel_diamond((5,5)),
                                 dist_threshold=4, min_overlap=10,
                                grayscale_threshold=0.3)
        adapt_bdry = falses(size(preprocessed_image))
        for tile in filtered_tiles
            adapt_bdry[tile...] = isboundary((labels_map(adapt_watershed) .== 0)[tile...])
        end
    end

    # To do: 
    # - determine best grayscale enhancement factors for brightening
    @info "Segmentation step 3: Enhanced grayscale k-means"
    begin
        segmentation_intersection = closing(kmeans_binarized, strel_diamond((3, 3))) .* adaptive_thresh
        boundary_intersection = (kmeans_bdry  .> 0 ) .&& (adapt_bdry .> 0)
        
        morphed_grayscale = grayscale_reconstruction(
            preprocessed_image,
            segmentation_intersection,
            boundary_intersection,
            land_mask;
            brightness_threshold = s.grayscale_floe_threshold)

        # New kmeans binarization, accounting for boundaries in initial segmentation
        kmeans_refined =
            kmeans_binarization(
                morphed_grayscale,
                falsecolor_image;
                k = 3,
                cluster_selection_algorithm=s.bright_ice_detector)
        kmeans_refined[boundary_intersection] .= 0

        kmeans_refined .= clean_binary_floes_mask_fill(
            kmeans_refined, cloud_mask, prelim_ice_mask; 
            max_fill=s.maximum_fill_size, min_floe_size=s.minimum_floe_size, opening_strel=strel_diamond((3,3)))

    end

    @info "Segmentation step 4: Morphological floe splitting"
    begin
        separated_floes = morph_split_floes(label_components(kmeans_refined), max_depth=5)
    end

    # @info "Image gradient filtering"

    # Original algorithm does k-means again after zero-ing out the water pixels
    # Then watershed -> remove smaller than 3 km2 -> hbridge -> opening with 5,5 strel -> fill holes.
    # FS pipeline: Segmentation F
    # BG pipeline: Regularize + Get Final
    # Side by side comparison -- which does a better job, given a good initial binarization? 

    # Final step in this version: select low-solidity objects, and break the bridges
    # Need some confidence that we are separating out the filaments to throw away
    return SegmentedImage(truecolor_image, separated_floes)
end


###### Tracking #######
# Specify defaults for the various filter functions


###### Helper Functions ######
# Functions needed for specific parts of the workflow. 
# TODO: ice_water_mask should be an IceDetectionAlgorithm and could be in the main code base
# TODO: ice_water_mask could also take a masked TC image, so we don't have to do the apply! functions again.
"""
    ice_water_mask(tc_img_, cloud_mask, land_mask; b1_ice_min = 75/255)

Identify sea ice pixels using the Band 1 reflectance.
""" 
function ice_water_mask(truecolor_image, cloud_mask, land_mask; b1_ice_min = 75/255) 
    tc_img = RGB.(truecolor_image)
    apply_cloudmask!(tc_img, cloud_mask)    
    apply_landmask!(tc_img, land_mask .> 0)
    
    banddata = red.(tc_img)
    
    ##### replace with binarize function when written ######
    edges, bincounts = build_histogram(banddata, 64; minval=0, maxval=1)
    ice_peak = get_ice_peaks(edges, bincounts;
        possible_ice_threshold=b1_ice_min,
        minimum_prominence=0.01,
        window_size=3)
    
    thresh = 0.5 * (b1_ice_min + ice_peak)
    return banddata .> thresh
end

####### Morphological cleanup ########
"""
    clean_binary_floes(bw_img; min_opening_area=50, min_object_size=16)
    clean_binary_floes(bw_img, tiles; min_opening_area=50, min_object_size=16)

    Use morphological operations to separate connect floes, fill holes, and remove small speckles.
    - `bw_img`: Binary mask with objects of interest = 1 or true
    - `tiles`: (optional) TileIterator, e.g. [`get_tiles`](@ref)
    - `min_opening_area`: `min_area` parameter sent to the `area_opening` function from ImageMorphology. Used 
      as in input to `imfill` as well.
    - `min_object_area`: Minimum size object to retain in `bw_img`. 

"""
function clean_binary_floes(
    bw_img::AbstractArray{Bool};
    min_opening_area::Int=50,
    min_object_size::Int=16
)
    img_opened = area_opening(bw_img; min_area=min_opening_area) |> hbreak
    img_filled = branch(img_opened) |> bridge
    img_filled = .!imfill(.!img_filled, (0, min_opening_area))
    diff_matrix = img_opened .!= img_filled
    cleaned_img = bw_img .|| diff_matrix
    return imfill(cleaned_img, (0, min_object_size))
end

function clean_binary_floes(
    bw_img::AbstractArray{Bool},
    tiles;
    min_opening_area::Int=50,
    min_object_size::Int=16
)
    out = falses(size(bw_img))
    for tile in tiles
        out[tile...] .= clean_binary_floes(bw_img[tile...];
                            min_opening_area=min_opening_area,
                            min_object_size=min_object_size)
    end
    return out
end

##### Watershed transformation #######
"""
   watershed_transform(binary_img, img; strel=strel_diamond((3,3)), dist_threshold=4)
   watershed_transform(binary_img, img, tiles; strel=strel_diamond((3,3)), dist_threshold=4)

    Carry out a watershed transform of `binary_img` by computing the distance transform,
    selecting marker regions with distance to background greater than `dist_threshold`, then eroding
    the marker regions using `strel`. Following the watershed transform, the background of `binary_img`
    is imposed and a SegmentedImage is generated using `img`. Optionally, supply a TiledIterator and 
    carry out the transform only on the listed tiles.

"""
function watershed_transform(
    binary_floes,
    img;
    strel=strel_diamond((5,5)),
    dist_threshold=4
)
    bw = .!binary_floes
    dist = .- distance_transform(feature_transform(bw))
    markers = erode(dist .< -1 * dist_threshold, strel) |> label_components
    labels = labels_map(watershed(dist, markers))
    labels[bw] .= 0

    return SegmentedImage(img, labels)
end

function watershed_transform(
    binary_floes,
    img,
    tiles;
    strel=strel_diamond((5,5)),
    dist_threshold=4,
    min_overlap=10,
    grayscale_threshold=0.1
)
    labels = zeros(Int64, size(binary_floes))
    for tile in tiles
        wseg = watershed_transform(binary_floes[tile...], img[tile...]; strel=strel, dist_threshold=dist_threshold)
        labels[tile...] = labels_map(wseg)
    end
    wseg = SegmentedImage(img, labels)
    stitched_labels = stitch_clusters(wseg, tiles, min_overlap, grayscale_threshold) # TODO: Add method to allow labels map for stitch clusters
    stitched_labels[.!binary_floes] .= 0 # TODO: Check to see if the background needs to be re-masked

    return SegmentedImage(img, stitched_labels)
end

# Component properties can go into the regionprops file
"""
    component_boundary_mean(indexmap, img, strel)

Compute the mean of `img` within the external boundary of objects
in `indexmap`. We define the external boundary of I as the set
`dilate(I) ∖ I`.
"""
function component_boundary_mean(indexmap, img, strel)
    labels = unique(indexmap)
    results = Dict()
    for l in labels
        idx = indexmap .== l
        bdry = dilate(idx, strel) .- idx
        results[l] = mean(vec(img[bdry .> 0]))
    end
    return results
end 

# We can add this to the regionprops table if wanted

"""
    component_solidities(indexmap)

Compute the solidity (area/convex area) for each labeled region in indexmap. 

"""
function component_solidities(indexmap)
    ca = component_convex_areas(indexmap)
    a = component_lengths(indexmap)
    s = Dict(r => (r == 0 ? 0 : a[r] / ca[r]) for r in keys(ca))
    return s
end

"""
    component_circularities(indexmap)

Compute the circularity (4πA/P^2)for each labeled region in indexmap. 

"""
function component_circularities(indexmap)
    p = component_perimeters(indexmap)
    a = component_lengths(indexmap)
    c = Dict(r => (r == 0 ? 0 : 4*π*a[r] / p[r]^2) for r in keys(a))
    return c
end

"""
    grayscale_reconstruction(preprocessed_image, binarized_image, boundary_intersection;
            brightness_threshold = s.grayscale_floe_threshold,
            min_area_opening=20,
            strel_dilation = strel_diamond((5, 5))
            )

    Use the information in a binarized image and detected boundaries to enhance a grayscale image.
    The image is processed by darkening leads, then using reconstruction to smooth regions. This method
    is based on the LopezAcosta2021 Fram Strait pipeline.
"""
function grayscale_reconstruction(preprocessed_image, binarized_image, boundary_intersection, land_mask;
    brightness_threshold = 0.4,
    brightening_factor = 0.3,
    min_area_opening=20,
    strel_dilation = strel_diamond((5, 5))
    )

    # Enhance grayscale image - can be a function
    brightened_grayscale = deepcopy(preprocessed_image) 
    brightened_grayscale[brightened_grayscale .< brightness_threshold] .= 0
    brightened_grayscale .= adjust_histogram(
        brightened_grayscale .* brightening_factor .+ preprocessed_image, LinearStretching()
        )
    apply_landmask!(brightened_grayscale, land_mask)
                                     # could use a robust linear stretching with percentiles



    ice_leads = .!boundary_intersection .* binarized_image
    ice_leads .= .!area_opening(ice_leads; min_area=min_area_opening, connectivity=2)


    brightened_dilated = dilate(brightened_grayscale, strel_dilation)

    mreconstruct!(
        dilate, brightened_dilated,
        complement.(brightened_dilated),
        complement.(brightened_grayscale)
    )
    return adjust_histogram(brightened_dilated .* ice_leads, LinearStretching())
end

"""
    clean_binary_floes_mask_fill

Fill small holes in the binary image if the hole represents clouds or ice. Use morphological operations (hbreak, opening)
to clean up the binary image for floe detection.
"""
function clean_binary_floes_mask_fill(
    binary_image, cloud_mask, prelim_ice_mask; 
    max_fill=64, min_floe_size=100, opening_strel=strel_diamond((3,3))
)
    binary_image .= imfill(
                        hbreak(binary_image),
                        (0, min_floe_size)
                    )
    # Fill moderately small holes if they are clouds or ice
    small_holes = .!imfill(.!binary_image, (0, max_fill)) .&& .! binary_image
    # To do: Find small holes with no neighbors
    
    small_holes .= small_holes .&& (cloud_mask .|| prelim_ice_mask)
    binary_image[small_holes .> 0] .= 1    
    binary_image .= opening(binary_image, opening_strel)
    binary_image .= imfill(binary_image, (0, 100)) 
    return binary_image
end

using DataFrames
import IceFloeTracker.Morphology: _generate_se!
function se_disk(r)
    se = [sum(c.I) <= r for c in CartesianIndices((2*r + 1, 2*r + 1))]
    _generate_se!(se)
    return se
end

"""
    morph_split_floes(indexmap; max_depth=5, filter_function)

Split floes using morphological opening. For each connected component in `binary_img`,
try applying morphological opening up to a maximum depth of `max_depth` until the object
is split into more than one piece. Optionally, supply a list of labels and only apply the
splitting to those labels.

Filter function can use the LinearThresholdFunctions from the tracking library. For now, filter function
needs to be a function that takes circularity and solidity as inputs, and produces a list of labels which 
are too low to count as floes.

Options for later:
- inplace version
- function for choosing which labels to split
- option to keep labels intact
"""
function morph_split_floes(labeled_array;
        max_depth=5,
        min_area=100,
        filter_function=(c, s) -> (c < 0.6) | (s < 0.85),
        max_iter=10)

    # relabel components to remove gaps in numbering
    indexmap = label_components(labeled_array)
    out = deepcopy(indexmap)
    
    label_offset = maximum(indexmap)
    boxes = component_boxes(indexmap)
    indices = component_indices(indexmap)
    masks = component_floes(indexmap)

    # apply opening with radius r = 1:max_depth until either the
    # shape is broken into parts or you reach max depth. If max depth
    # is reached and nothing changes, return the original mask.
    function morph_split(mask, max_depth)
        for r in 1:max_depth
            # d = 2*r + 1 # This would be if we wanted instead to use the strel_diamond or box
            update_mask = opening(mask, se_disk(r)) |> label_components
            maximum(update_mask) > 1 && return update_mask .> 0
        end
        return mask
    end

    # TODO: Make this an input function, not hard coded
    areas = component_lengths(indexmap)
    perimeters = component_perimeters(indexmap)
    convex_areas = component_convex_areas(indexmap)
    labels = filter(r -> r != 0, intersect(keys(areas), keys(perimeters), keys(convex_areas)))
    
    circularities = Dict(r => 4 * pi * areas[r] / perimeters[r]^2 for r in labels if r != 0)
    solidities = Dict(r => areas[r] / convex_areas[r] for r in labels if r != 0)

    split_labels = [r for r in labels if filter_function(circularities[r], solidities[r])]

    n_regions = length(unique(out))
    n_updated = 0
    split_labels_new = []

    count = 0
    
    while (n_regions != n_updated) && (count < max_iter)
        n_regions = length(unique(out))
        for r in split_labels
            update_floe = morph_split(masks[r], max_depth)
            masks[r] != update_floe && 
                begin
                    update_labels = label_components(update_floe)
                    update_areas = component_lengths(update_labels)
                    for rnew in keys(update_areas)
                        # imfill can work, but won't catch if a small object neighbors another.                
                        update_areas[rnew] < min_area && (update_labels[update_labels .== rnew] .= 0)
                    end
                     # re-label for dropped objects
                    update_labels = label_components(update_floe)
                    update_labels[update_labels .> 0] .+= label_offset
                    # update offset based on the new labels
                    label_offset += length(unique(update_labels))
    
                    # add new shape(s) to the indexmap
                    out[indices[r]] .= 0 # Remove original floe
                    out[boxes[r]] .+= update_labels # New floe has at most equal area to the original, so there shouldn't be overlap
    
                    # check if new components need to be added to the update list
                    # TODO: make this all part of the filter function
                    a = component_lengths(update_labels)
                    p = component_perimeters(update_labels)
                    ca = component_convex_areas(update_labels)
                    _labels = [r for r in unique(update_labels) if (r != 0) && (r ∈ keys(a)) && (r ∈ keys(p)) && (r ∈ keys(ca))]
                    for rnew in _labels 
                        (rnew != 0) && begin
                            c = 4 * pi * a[rnew] / p[rnew]^2
                            s = a[rnew] / ca[rnew]
                            filter_function(c, s) && (push!(split_labels_new, rnew))
                        end
                    end
                end
        end
        # update the split_labels list and zero out the split_labels_new list
        count += 1
        n_updated = length(unique(out))

        # Update maps
        boxes = component_boxes(out)
        indices = component_indices(out)
        masks = component_floes(out)
        split_labels = [r for r in split_labels_new if r in keys(masks)]
        split_labels_new = []

        # print("Count: "*string(count)*", N labels "*string(n_regions)*", N updated "*string(n_updated)*"\n")
    end
    return label_components(out)
end

end
