module Watkins2026
import Images
using TiledIteration: TileIterator

import ..Filtering: nonlinear_diffusion, PeronaMalikDiffusion, unsharp_mask, channelwise_adapthisteq
import ..Morphology: hbreak, hbreak!, branch, bridge, fill_holes
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
    IceDetectionBrightnessPeaksMODIS721,
    stitch_clusters,
    regionprops
end

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
    coastal_buffer_structuring_element = strel_box((51,51))
    cloud_mask_algorithm = Watkins2025CloudMask()
    diffusion_algorithm = PeronaMalikDiffusion(lambda=0.1, kappa=0.1, niters=5, g="exponential")
    adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.99) # rblocks/cblocks not used yet -- add with CLAHE.jl
    unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
end

# TODO: Decide how to pass landmasks and cloudmasks around
function (p::Preprocess)(
    truecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}}, 
    falsecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},  
    landmask_image,
    tiles::TileIterator
)::AbstractArray{AbstractGray}
    
    proc_img = deepcopy(Gray.(truecolor_image))
    cloud_mask = p.cloud_mask_algorithm(falsecolor_image)
    
    _lm_temp = create_landmask(landmask_image, p.coastal_buffer_structuring_element)
    land_mask = _lm_temp.non_dilated

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
    begin
        proc_img = Gray.(IceFloeTracker.skimage.sk_exposure.equalize_adapthist(
            IceFloeTracker.ImageUtils.to_uint8(Float64.(proc_img) .* 255);
            # Using default: image size divided by 8. Update when CLAHE.jl available
            clip_limit = 1 - p.adapthisteq_params.clip,  # Equivalent to MATLAB's 'ClipLimit'
            nbins=p.adapthisteq_params.nbins,      # Number of histogram bins. 255 is used to match the default in MATLAB script
        ))
    end

    # Reapply masks, since sharpening and eq adjustment can bleed into neighboring pixels
    apply_cloudmask!(proc_img, cloud_mask)
    apply_landmask!(proc_img, land_mask)
    
    return proc_img
end

###### Segmentation ########


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
    ice_peak = IceFloeTracker.get_ice_peaks(edges, bincounts;
        possible_ice_threshold=b1_ice_min,
        minimum_prominence=0.01,
        window=3)
    
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
    tiles::TileIterator;
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
function watershed_transform(binary_floes, img; strel=strel_diamond((3,3)), dist_threshold=4)
    bw = .!binary_floes
    dist = .- distance_transform(feature_transform(bw))
    markers = erode(dist .< -1 * dist_threshold, strel) |> label_components
    labels = labels_map(watershed(dist, markers))
    labels[bw] .= 0
    return SegmentedImage(img, labels)
end

function watershed_transform(binary_floes, img, tiles; strel=strel_diamond((3,3)), dist_threshold=4, min_overlap=10)
    labels = zeros(Int64, size(binary_floes))
    for tile in tiles
        wseg = watershed_transform(binary_floes[tile...], img[tile...]; strel=strel, dist_threshold=dist_treshold)
        labels[tile...] = labels_map(wseg)
    end
    wseg = SegmentedImage(img, labels)
    stitched_labels = stitch_clusters(wseg, tiles, min_overlap) # TODO: Add method to allow labels map for stitch clusters
    stitched_labels[.!binary_floes] .= 0 # TODO: Check to see if the background needs to be re-masked
    return SegmentedImage(img, stitched_labels)
end

# Component properties
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