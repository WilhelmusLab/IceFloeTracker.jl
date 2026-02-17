module Watkins2026
using Images # Find individual functions to import later

import ..Filtering: nonlinear_diffusion, PeronaMalikDiffusion, unsharp_mask, channelwise_adapthisteq
import ..Morphology: hbreak, hbreak!, branch, bridge, _generate_se!
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
    prelim_ice_threshold = 0.3
    bright_ice_detector = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01
    )
    grayscale_floe_threshold = 0.4 # Brighten pixels brighter than this threshold 
    adaptive_binarization_settings = (minimum_window_size=400, threshold_percentage=0, minimum_brightness=100/255)
    watershed_strel = se_disk(5) 
    # cleanup parameters: can be named tuple, after writing self-contained cleanup function
    minimum_floe_size = 100
    maximum_floe_size = 50000 # how big should we allow?
    maximum_fill_size = 64
    morph_max_depth = 10
    #kmeans_binarization_settings = (k=4, )
end

function (s::Segment)(truecolor_image, falsecolor_image, land_mask_image)

    n, m = size(truecolor_image)
    tile_size_pixels = s.tile_size_pixels
    tile_size_pixels > maximum([n, m]) && begin
        @warn "Tile size too large, defaulting to image size"
        tile_size_pixels = minimum([n, m])
    end
    tiles = get_tiles(truecolor_image, tile_size_pixels)

    cloud_mask = s.cloud_mask_algorithm(falsecolor_image)
    
    _lm_temp = create_landmask(land_mask_image, s.coastal_buffer_structuring_element)
    land_mask = _lm_temp.non_dilated
    coastal_buffer = _lm_temp.dilated
    
    # Ice-water mask from the red channel reflectance
    # TODO: Add as an IceDetectionAlgorithm to use as an import
    prelim_ice_mask = ice_water_mask(truecolor_image, cloud_mask, coastal_buffer, tiles; b1_ice_min=s.prelim_ice_threshold)

    # Enhanced grayscale image for segmentation
    preprocessed_image = s.preprocessing_function(truecolor_image, land_mask, cloud_mask, tiles)

    # Select a subset of the tiles to preprocess by checking for clear-sky ice pixels outside the 
    # coast buffer
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > s.min_tile_ice_pixel_count,
        tiles);

    @info "Segmentation method 1: First K-means workflow"
    # TODO: Test impact of the three different "floe brightening" methods
    # - Morphological residue method from LopezAcosta2019Tiling
    # - Discriminate-ice-water from LopezAcosta2019
    # - Discriminate-ice-water from the LASW2019 paper
    # - Gamma threshold method from LopezAcosta2019
    # - Make the kmeans_binarized a binarization algorithm
    begin
        fc_masked = apply_cloudmask(RGB.(falsecolor_image), cloud_mask) # Error in cloudmask function with RGBA
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
                ) .> 0 # returns Gray, so we convert it here
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
    @info "Preprocess 2: Enhanced grayscale"
    begin
        # Simple join of segmentation results: closing (dilate then erode) of the k-means results
        # then multiplication by the prelim ice mask and the adaptive threshold mask
        # TODO: Segmentation intersection could be a function.
        segmentation_intersection = closing(kmeans_binarized, strel_diamond((3, 3))) .* adaptive_thresh .* prelim_ice_mask
        boundary_intersection = (kmeans_bdry  .> 0 ) .&& (adapt_bdry .> 0)
        
        # The grayscale reconstruction method comes from the LopezAcosta2019 code
        # The idea is that we use the segmentation results from the initial processing to reconstruct
        # the image, smearing out image details based on where the mask is. It also uses a 
        # threshold to get an extra segmentation method for brightening the image.
        morphed_grayscale = grayscale_reconstruction(
            preprocessed_image,
            segmentation_intersection,
            boundary_intersection,
            land_mask;
            brightness_threshold = s.grayscale_floe_threshold)
        morphed_grayscale[boundary_intersection] .= 1 # Testing applying this before instead of after the k-means step
    end

    @info "Segmentation 3: K-means of morphed grayscale image"
    begin
        # Then, we use k-means on the image again.
        # New kmeans binarization, accounting for boundaries in initial segmentation
        kmeans_refined =
            kmeans_binarization(
                morphed_grayscale,
                fc_masked,
                filtered_tiles;
                k = 3,
                cluster_selection_algorithm=s.bright_ice_detector)
        # kmeans_refined[boundary_intersection] .= 0 # Q: better before or after?

        # The third application of the morphological cleanup removes small floes
        # and fills holes with certain criteria. Needs to be updated
        kmeans_refined .= clean_binary_floes_mask_fill(
            kmeans_refined, cloud_mask, prelim_ice_mask; 
            max_fill=s.maximum_fill_size, min_floe_size=s.minimum_floe_size, opening_strel=strel_diamond((3,3)))
        
        # Add special treatment for floes that are in the cloud mask bounday. For these floes, 
        # we can take the intersection between the convex 

        apply_landmask!(kmeans_refined, land_mask)
        apply_cloudmask!(kmeans_refined, cloud_mask)
    end

    return SegmentedImage(truecolor_image, label_components(kmeans_refined)) # Temp early exit
    @info "Floe splitting"
    # TODO: Floe splitting can be an input to the function.
    # TODO: Add parameters to main function call
    begin
        # morph split goes from binarized to index map. Update the outputs to be consistent across floe split functions.
        separated_floes = morph_split_floes(label_components(kmeans_refined))
        
        # separated_floes = label_components(kmeans_refined) # temp test if morph split helps
        # watershed split floes goes from indexmap to segmented image
        watershed_floes = watershed_split_floes(separated_floes, truecolor_image;
                            strel=s.watershed_strel, min_circularity=0.6, min_solidity=0.85)
    end

    @info "Final cleanup"

    # Enforce the min floe size, potential to remove floe candidates if other criteria aren't met
    indexmap = labels_map(watershed_floes)
    indexmap .= indexmap .* opening(indexmap .> 0, se_disk(2)) # size of the disk can be a parameter
    areas = component_lengths(indexmap)
    indices = component_indices(indexmap)
    final_indexmap = zeros(Int64, size(indexmap))

    # could set up a function that collects the labels of the 
    # floes that fulfil the final criteria
    # - aspect ratios, circularity, solidity, gradient ratio

    label = 1
    for r in keys(areas)
        (r != 0) && (
            s.minimum_floe_size <= areas[r] <= s.maximum_floe_size) && begin
                final_indexmap[indices[r]] .= label
                label += 1
        end 
    end
    
    # Check brightness, edge contrast first
    return SegmentedImage(truecolor_image, final_indexmap)
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

function ice_water_mask(truecolor_image, cloud_mask, land_mask, tiles; b1_ice_min = 75/255) 
    out = falses(size(truecolor_image))
    for tile in tiles
        out[tile...] .= ice_water_mask(truecolor_image[tile...], cloud_mask[tile...], land_mask[tile...]; b1_ice_min=b1_ice_min)
    end
    return out
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

""" # TODO: Make the "marker_selection_function" an input, with the idea that it's a function of the distance transform
function watershed_transform(
    binary_floes,
    img;
    strel=strel_diamond((5,5)),
    dist_threshold=4
)
    bw = .!binary_floes
    dist = .- distance_transform(feature_transform(bw))
    markers = erode(dist .< -dist_threshold, strel) |> label_components
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
    labels = intersect(keys(a), keys(ca))
    s = Dict(r => (r == 0 ? 0 : a[r] / ca[r]) for r in labels)
    return s
end

"""
    component_circularities(indexmap)

Compute the circularity (4πA/P^2)for each labeled region in indexmap. 

"""
function component_circularities(indexmap)
    p = component_perimeters(indexmap)
    a = component_lengths(indexmap)
    labels = intersect(keys(a), keys(p))
    c = Dict(r => (r == 0 ? 0 : 4*π*a[r] / p[r]^2) for r in labels)
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
""" # TODO: Set up function to select floes that intersect with cloud mask, and test the hole and gap-filling methods.
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
    # Seems to me like this could be done with the distance transform or with morphology
    small_holes .= small_holes .&& (cloud_mask .|| prelim_ice_mask)
    binary_image[small_holes .> 0] .= 1    
    binary_image .= opening(binary_image, opening_strel)
    binary_image .= imfill(binary_image, (0, 100)) 
    return binary_image
end

"""
    se_disk(r)

    Generate an approximately circular structuring element with radius r. For small r, this will be somewhat diamond-shaped.

""" # #TODO add simple example to docs, add to special strels, and in future, optimize the extreme filter for this shape
function se_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    return sqrt.(se) .<= r
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
""" # TODO: Add this function to a new floe-splitting segmentationl.jl file, and add some tests for cases (like linear features)
function morph_split_floes(labeled_array;
        max_depth=20,
        min_area=100,
        min_circularity=0.6,
        min_solidity=0.85,
        max_iter=10)

    # relabel components to remove gaps in numbering
    indexmap = label_components(labeled_array)
    out = deepcopy(indexmap)
    boxes = component_boxes(indexmap)
    indices = component_indices(indexmap)
    masks = component_floes(indexmap)

    # Split the mask by opening with a circular SE until the floes separate.
    # If it isn't split after reaching the effective radius, quit.
    function morph_split(mask, max_depth)
        _max_depth = round(Int, sqrt(sum(mask)/pi))
        # check for oblong shapes
        n, m = size(mask)
        _max_depth = minimum([n-1, m-1, _max_depth])
        for r in 1:_max_depth
            new_labels = opening(mask, se_disk(r)) |> label_components
            length(unique(new_labels)) > 2 && return new_labels .> 0
        end
        return nothing
    end

    # Get the list of candidates for splitting by checking area, circularity, and solidity.
    # Evaluates solidity lazily since that operation is slow.
    function get_candidate_labels(indexmap, min_area, min_circularity, min_solidity)
        areas = component_lengths(indexmap)
        perimeters = component_perimeters(indexmap)
        masks = component_floes(indexmap)
        # convex_areas = component_convex_areas(indexmap)
        labels = filter(r -> r != 0, intersect(keys(areas), keys(perimeters)))
        
        circularities = Dict(r => 4 * pi * areas[r] / perimeters[r]^2 for r in labels if r != 0)
        
        # Only evaluate convex area if the area and circularity pass the thresholds for splitting
        split_labels = Int64[]
        for r in labels
            (r != 0) && (areas[r] > 1.5 * min_area) && (circularities[r] < min_circularity) && begin
                ca = component_convex_areas(Int64.(masks[r]))
                s = areas[r] / ca[1]
                s < min_solidity && push!(split_labels, r)
            end
        end
        return split_labels
    end

    n_regions = length(unique(out))
    n_updated = 0
    split_labels = get_candidate_labels(indexmap, min_area, min_circularity, min_solidity)
    count = 0
    while (n_regions != n_updated) && (count < max_iter)
        n_regions = length(unique(out))
        for r in split_labels
            update_floe = morph_split(masks[r], max_depth)
            !isnothing(update_floe) && begin
                update_labels = label_components(update_floe)
                update_areas = component_lengths(update_labels)
                
                for rnew in keys(update_areas)
                    # imfill works for fully separated components, but will miss if components are touching.                
                    update_areas[rnew] < min_area && (update_labels[update_labels .== rnew] .= 0)
                end

                # Remove original floe and overwrite with new floe
                out[indices[r]] .= 0 
                out[boxes[r]] .+= update_labels 
            end
        end
        
        count += 1
        out .= label_components(out)
        n_updated = length(unique(out))
        n_regions != n_updated && begin 
            # Select candidate labels which were successfully split in the morph split step.
            split_labels = get_candidate_labels(out, min_area, min_circularity, min_solidity)
            ### There's an error here that doesn't make sense to me, since this code worked in a standalone script.
            ### If I can make it work, it should save at least some time for processing.
            # split_labels = [r for r in get_candidate_labels(out, min_area, min_circularity, min_solidity) if length(unique(out[indices[r]])) > 1]
            
            # One thing that slows this down is having to re-write the array all the time. It would be faster to add entries to the boxes
            # rather than start afresh each time.
            boxes = component_boxes(out)
            indices = component_indices(out)
            masks = component_floes(out)
        end
    end
    return label_components(out)
end

"""
    watershed_split_floes(indexmap, img; strel=strel_box((5,5)))

Attempt to separate floes using a watershed transformation. Finds markers using an adaptive method,
taking each reach and checking each depth from 1 to max_depth until the floe is separated at least once.
Carry out a watershed transform. If floes are separated, then overwrite the original indexmap. Otherwise,
don't change the floe.
"""
function watershed_split_floes(indexmap, img; strel=strel_box((5,5)), min_circularity=0.6, min_solidity=0.5)

    function dist_split(distances, max_depth)
        # TODO Make it possible to label and keep the distinct components
        # This could look like saving it each time the component separates
        # Or it could look like the E-P method, where we mark when a label disappears.
        for r in 1:max_depth
            update_mask = (distances .< -r) |> label_components
            maximum(update_mask) > 1 && return distances .< -r
        end
        return zeros(size(distances))
    end

    C = component_circularities(indexmap)
    S = component_solidities(indexmap)
    indices = component_indices(indexmap)
    masks = component_floes(indexmap)
    boxes = component_boxes(indexmap)
    
    labels = intersect(keys(C), keys(S), keys(indices))

    adjust_labels = [r for r in labels if (C[r] < min_circularity || S[r] < min_solidity)]

    # Binary image is the inverted map
    binary_img = indexmap .== 0

    # Initialize markers with the nonzero regions in the indexmap
    init_dist = .- distance_transform(feature_transform(binary_img))
    adjust_floes = zeros(Int64, size(indexmap));
    markers = zeros(Int64, size(indexmap));
    
    # Only transform floes that fail the circ/solidity check
    for r in adjust_labels
        adjust_floes[indices[r]] .= (indexmap .> 0)[indices[r]]
        update_mask = dist_split(init_dist[boxes[r]] .* masks[r], 5)
        markers[boxes[r]] .+= update_mask
    end

    bw_img = complement.(mreconstruct(dilate, markers, adjust_floes)) .> 0;
    updated_dist = .- distance_transform(feature_transform(bw_img));
    new_labels = labels_map(watershed(updated_dist, label_components(markers)));
    new_labels[bw_img .> 0] .= 0
    
    # Find where the new boundaries separate regions and exclude where they mark the boudnary with background
    watershed_bdry = (isboundary(new_labels) .> 0) .&& (isboundary(new_labels .> 0) .> 0)
    
    updated_indexmap = deepcopy(indexmap)
    updated_indexmap[new_labels .> 0] .= new_labels[new_labels .> 0]
    updated_indexmap[watershed_bdry] .= 0

    # Option later to use different marker selection method, perhaps using image gradient
    # Note that the watershed transform takes a binary image as input
    # watershed_segments = watershed_transform(watershed_floes .> 0, img; strel=strel)
    # watershed_indexmap = labels_map(watershed_segments)
    # watershed_bdry = isboundary(watershed_indexmap)
    # watershed_indexmap[watershed_bdry .> 0] .= 0

    # If the watershed transform separated the object into multiple parts, overwrite the original matrix
    # Using the component indices makes this much faster -- you only have to search the matrix for locations
    # one time.
    # for r in adjust_labels
    #     length(unique(watershed_indexmap[indices[r]])) > 1 && begin
    #         indexmap[indices[r]] .= watershed_indexmap[indices[r]]
    #     end
    # end

    return SegmentedImage(img, updated_indexmap)
#### Notes for next steps
# - Fix the adaptive marker finder. It needs to get at least a few levels.
# - Fix the issue where we get a line inbetween floes. Why is that happening?
# - Fix the issue where large objects are not being removed.
# - Use tracking between the pairs of images to find matched floes
# - Use image gradients to filter cloud / landfast segments
# - Seems like the largest shapes didn't get removed like they should have
# - Use mismatched near-time floes to refine (e.g., average images after optimal alignment)
# - Use best matches across pairs to estimate location of other pairs
end

end
