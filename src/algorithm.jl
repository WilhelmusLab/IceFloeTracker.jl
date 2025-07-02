abstract type IceFloeSegmentationAlgorithm end

struct LopezAcosta2019 <: IceFloeSegmentationAlgorithm
    landmask_structuring_element::AbstractMatrix{Bool}
end

function LopezAcosta2019(; landmask_structuring_element=make_landmask_se())
    return LopezAcosta2019(landmask_structuring_element)
end

function (p::LopezAcosta2019)(
    truecolor_image::T, falsecolor_image::T, landmask_image::U
) where {T<:Matrix{RGB{Float64}},U<:AbstractMatrix}
    @info "building landmask"
    landmask_imgs = create_landmask(landmask_image, p.landmask_structuring_element)

    @info "Building cloudmask"
    cloudmask = create_cloudmask(falsecolor_image)

    # 2. Intermediate images
    @info "Finding ice labels"
    ice_labels = IceFloeTracker.find_ice_labels(falsecolor_image, landmask_imgs.dilated)

    @info "Sharpening truecolor image"
    # a. apply imsharpen to truecolor image using non-dilated landmask
    sharpened_truecolor_image = IceFloeTracker.imsharpen(
        truecolor_image, landmask_imgs.non_dilated
    )
    # b. apply imsharpen to sharpened truecolor img using dilated landmask
    sharpened_gray_truecolor_image = IceFloeTracker.imsharpen_gray(
        sharpened_truecolor_image, landmask_imgs.dilated
    )

    @info "Normalizing truecolor image"
    normalized_image = IceFloeTracker.normalize_image(
        sharpened_truecolor_image, sharpened_gray_truecolor_image, landmask_imgs.dilated
    )

    # Discriminate ice/water
    @info "Discriminating ice/water"
    ice_water_discrim = IceFloeTracker.discriminate_ice_water(
        falsecolor_image, normalized_image, copy(landmask_imgs.dilated), cloudmask
    )

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    segA = IceFloeTracker.segmentation_A(
        IceFloeTracker.segmented_ice_cloudmasking(ice_water_discrim, cloudmask, ice_labels)
    )

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = IceFloeTracker.segmentation_B(sharpened_gray_truecolor_image, cloudmask, segA)

    # Process watershed in parallel using Folds
    @info "Building watersheds"
    # container_for_watersheds = [landmask_imgs.non_dilated, similar(landmask_imgs.non_dilated)]

    watersheds_segB = [
        IceFloeTracker.watershed_ice_floes(segB.not_ice_bit),
        IceFloeTracker.watershed_ice_floes(segB.ice_intersect),
    ]
    watersheds_segB_product = IceFloeTracker.watershed_product(watersheds_segB...)

    # segmentation_F
    @info "Segmenting floes part 3/3"
    segF = IceFloeTracker.segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB_product,
        ice_labels,
        cloudmask,
        landmask_imgs.dilated,
    )

    @info "Labeling floes"
    labeled_floes = label_components(segF)

    # TODO: return ImageSegmentation.jl-style results
    return labeled_floes
end

@kwdef struct LopezAcosta2019Tiling <: IceFloeSegmentationAlgorithm
    # Landmask parameters
    landmask_structuring_element::AbstractMatrix{Bool} = make_landmask_se()

    # Tiling parameters
    tile_rblocks::Integer = 8
    tile_cblocks::Integer = 8

    # Ice labels thresholds
    ice_labels_prelim_threshold::Float64 = 110.0
    ice_labels_band_7_threshold::Float64 = 200.0
    ice_labels_band_2_threshold::Float64 = 190.0
    ice_labels_ratio_lower::Float64 = 0.0
    ice_labels_ratio_upper::Float64 = 0.75

    # Adaptive histogram equalization parameters
    adapthisteq_white_threshold::Float64 = 25.5
    adapthisteq_entropy_threshold = 4
    adapthisteq_white_fraction_threshold::Float64 = 0.4

    # Gamma parameters
    gamma::Float64 = 1
    gamma_factor::Float64 = 1
    gamma_threshold::Float64 = 220

    # Unsharp mask parameters
    unsharp_mask_radius::Int = 10
    unsharp_mask_amount::Float64 = 2.0
    unsharp_mask_factor::Float64 = 255.0

    # Brighten parameters
    brighten_factor::Float64 = 0.1

    # Preliminary ice mask parameters
    prelim_icemask_radius::Int = 10
    prelim_icemask_amount::Int = 2
    prelim_icemask_factor::Float64 = 0.5

    # Main ice mask parameters
    icemask_band_7_threshold::Int = 5
    icemask_band_2_threshold::Int = 230
    icemask_band_1_threshold::Int = 240
    icemask_band_7_threshold_relaxed::Int = 10
    icemask_band_1_threshold_relaxed::Int = 190
    icemask_possible_ice_threshold::Int = 75
    icemask_n_clusters::Int = 3
end

function (p:LopezAcosta2019Tiling)(;
    truecolor_image::T, falsecolor_image::T, landmask_image::U
) where {T<:Matrix{RGB{Float64}},U<:AbstractMatrix}

    # Invert the landmasks – in the tiling version of the code, 
    # the landmask is expected to be the other polarity compared with
    # the non-tiling version.
    @info "building landmask"
    landmask_imgs = create_landmask(landmask_image, p.landmask_structuring_element)
    landmask = (dilated=.!landmask_imgs.dilated,)

    @info "Remove alpha channel if it exists"
    rgb_truecolor_img = RGB.(truecolor_image)
    rgb_falsecolor_img = RGB.(falsecolor_image)

    @info "Get tile coordinates"
    tiles = IceFloeTracker.get_tiles(
        rgb_truecolor_img; rblocks=p.tile_rblocks, cblocks=p.tile_cblocks
    )
    @debug tiles

    @info "Set ice labels thresholds"
    ice_labels_thresholds = (
        prelim_threshold=p.ice_labels_prelim_threshold,
        band_7_threshold=p.ice_labels_band_7_threshold,
        band_2_threshold=p.ice_labels_band_2_threshold,
        ratio_lower=p.ice_labels_ratio_lower,
        ratio_upper=p.ice_labels_ratio_upper,
        use_uint8=true,
    )
    @debug ice_labels_thresholds

    @info "Set adaptive histogram parameters"
    adapthisteq_params = (
        white_threshold=p.adapthisteq_white_threshold,
        entropy_threshold=p.adapthisteq_entropy_threshold,
        white_fraction_threshold=p.adapthisteq_white_fraction_threshold,
    )
    @debug adapthisteq_params

    @info "Set gamma parameters"
    adjust_gamma_params = (
        gamma=p.gamma, gamma_factor=p.gamma_factor, gamma_threshold=p.gamma_threshold
    )
    @debug adjust_gamma_params

    @info "Set structuring elements"
    # This isn't tunable in the underlying code right now, 
    # so just use the defaults
    structuring_elements = IceFloeTracker.structuring_elements
    @debug structuring_elements

    @info "Set unsharp mask params"
    unsharp_mask_params = (
        radius=p.unsharp_mask_radius,
        amount=p.unsharp_mask_amount,
        factor=p.unsharp_mask_factor,
    )
    @debug unsharp_mask_params

    @info "Set brighten factor"
    @debug p.brighten_factor

    @info "Set preliminary ice masks params"
    prelim_icemask_params = (
        radius=p.prelim_icemask_radius,
        amount=p.prelim_icemask_amount,
        factor=p.prelim_icemask_factor,
    )
    @debug prelim_icemask_params

    @info "Set ice masks params"
    ice_masks_params = (
        band_7_threshold=p.icemask_band_7_threshold,
        band_2_threshold=p.icemask_band_2_threshold,
        band_1_threshold=p.icemask_band_1_threshold,
        band_7_threshold_p.relaxed=icemask_band_7_threshold_relaxed,
        band_1_threshold_p.relaxed=icemask_band_1_threshold_relaxed,
        possible_ice_threshold=p.icemask_possible_ice_threshold,
        k=p.icemask_n_clusters, # number of clusters for kmeans segmentation
        factor=255, # normalization factor to convert images to uint8
    )
    @debug ice_masks_params

    @info "Segment floes"
    segmented_floes = IceFloeTracker.preprocess_tiling(
        n0f8.(rgb_falsecolor_img),
        n0f8.(rgb_truecolor_img),
        landmask,
        tiles,
        ice_labels_thresholds,
        adapthisteq_params,
        adjust_gamma_params,
        structuring_elements,
        unsharp_mask_params,
        ice_masks_params,
        prelim_icemask_params,
        p.brighten_factor,
    )

    @info "Label floes"
    labeled_floes = label_components(segmented_floes)

    # TODO: return ImageSegmentation.jl-style results
    return labeled_floes
end
