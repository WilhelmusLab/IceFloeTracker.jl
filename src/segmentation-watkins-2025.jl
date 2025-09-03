using IceFloeTracker: 
    LopezAcostaCloudMask,
    kmeans_segmentation

# This file implements an updated version of the lopez-acosta-2019-tiling segmentation algorithm
# following the IEEE TGRS calibration and validation paper. 
# I'm starting with John's adaptation of the tiling and non-tiling versions of the LA2019 algorithm.

abstract type IceFloeSegmentationAlgorithm end # defined elsewhere, can we remove it?

cloud_mask_thresholds = (
    prelim_threshold=53.0/255.,
    band_7_threshold=130.0/255.,
    band_2_threshold=169.0/255.,
    ratio_lower=0.0,
    ratio_offset=0.0,
    ratio_upper=0.53
)


# TBD
# Update to adapt the tile size to the image size. Tile size
@kwdef struct Watkins2025 <: IceFloeSegmentationAlgorithm
    tile_settings = (; rblocks=2, cblocks=2) 
    cloud_mask_thresholds = cloud_mask_thresholds
    adapthisteq_params = adapthisteq_params
    adjust_gamma_params = adjust_gamma_params
    structuring_elements = structuring_elements
    unsharp_mask_params = unsharp_mask_params
    ice_masks_params = ice_masks_params
    prelim_icemask_params = prelim_icemask_params
    brighten_factor = brighten_factor
end

# which things do we want to easily adapt between runs?
function Watkins2025(; landmask_structuring_element=make_landmask_se())
    return LopezAcosta2019(landmask_structuring_element) # Maybe move this to a setting
end

function (p::Watkins2025)(
    truecolor::T,
    falsecolor::T,
    landmask::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:AbstractMatrix}

    # Move these conversions down through the function as each step gets support for 
    # the full range of image formats
    truecolor_image = float64.(truecolor)
    falsecolor_image = float64.(falsecolor)
    landmask_image = float64.(landmask)

    @info "building landmask"
    landmask = create_landmask(landmask_image, p.landmask_structuring_element)

    @info "Building cloudmask"
    cm_algo = LopezAcostaCloudMask(p.cloud_mask_thresholds...)
    cloudmask = create_cloudmask(falsecolor_image, cm_algo)
    # TBD: Cloudmask cleanup. Fill small holes, remove speckle.

    @info "Tiling"
    tiles = get_tiles(truecolor; p.tile_settings...)
    # TBD: Set up tile function to pick the number of rblocks and cblocks that gets closest to having tiles 
    # of the specified size in pixels (100 km == 400 pixels)
    # TBD: set up a tile filter function to remove tiles with insufficient ocean pixels

    @info "Preprocessing"


    
    # 2. Intermediate images
    @info "Finding ice labels"
    ice_mask = find_ice_mask(falsecolor_image, landmask_imgs.dilated)

    @info "Sharpening truecolor image"
    # a. apply imsharpen to truecolor image using non-dilated landmask
    sharpened_truecolor_image = imsharpen(truecolor_image, landmask_imgs.non_dilated)
    # b. apply imsharpen to sharpened truecolor img using dilated landmask
    sharpened_gray_truecolor_image = imsharpen_gray(
        sharpened_truecolor_image, landmask_imgs.dilated
    )

    @info "Normalizing truecolor image"
    normalized_image = normalize_image(
        sharpened_truecolor_image, sharpened_gray_truecolor_image, landmask_imgs.dilated
    )

    # Discriminate ice/water
    @info "Discriminating ice/water"
    ice_water_discrim = discriminate_ice_water(
        falsecolor_image, normalized_image, copy(landmask_imgs.dilated), cloudmask
    )

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    segA = segmentation_A(
        segmented_ice_cloudmasking(ice_water_discrim, cloudmask, ice_mask)
    )

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = segmentation_B(sharpened_gray_truecolor_image, cloudmask, segA)

    # Process watershed in parallel using Folds
    @info "Building watersheds"
    # container_for_watersheds = [landmask_imgs.non_dilated, similar(landmask_imgs.non_dilated)]

    watersheds_segB = [
        watershed_ice_floes(segB.not_ice_bit), watershed_ice_floes(segB.ice_intersect)
    ]
    watersheds_segB_product = watershed_product(watersheds_segB...)

    # segmentation_F
    # TODO: @hollandjg find out why segF is more dilated
    @info "Segmenting floes part 3/3"
    segF = segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB_product,
        ice_mask,
        cloudmask,
        landmask_imgs.dilated,
    )

    @info "Labeling floes"
    labels = label_components(segF)

    # Return the original truecolor image, segmented
    segments = SegmentedImage(truecolor, labels)

    if !isnothing(intermediate_results_callback)
        segments_truecolor = SegmentedImage(truecolor, labels)
        segments_falsecolor = SegmentedImage(falsecolor, labels)
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            landmask_dilated=landmask_imgs.dilated,
            landmask_non_dilated=landmask_imgs.non_dilated,
            cloudmask=cloudmask,
            ice_mask=ice_mask,
            sharpened_truecolor_image=sharpened_truecolor_image,
            sharpened_gray_truecolor_image=sharpened_gray_truecolor_image,
            normalized_image=normalized_image,
            ice_water_discrim=ice_water_discrim,
            segA=segA,
            segB=segB,
            watersheds_segB_product=watersheds_segB_product,
            segF=segF,
            labels=labels,
            segment_mean_truecolor=map(
                i -> segment_mean(segments_truecolor, i), labels_map(segments_truecolor)
            ),
            segment_mean_falsecolor=map(
                i -> segment_mean(segments_falsecolor, i), labels_map(segments_falsecolor)
            ),
        )
    end
    return segments
end
