abstract type IceFloeTrackerAlgorithm end

struct LopezAcosta <: IceFloeTrackerAlgorithm
    landmask_structuring_element::AbstractMatrix{Bool}
end

function LopezAcosta(; landmask_structuring_element=make_landmask_se())
    return LopezAcosta(landmask_structuring_element)
end

function segment(
    truecolor::{<:AbstractArray{<:Union{AbstractRGB,TransparentRGB}}},
    falsecolor::{<:AbstractArray{<:Union{AbstractRGB,TransparentRGB}}},
    landmask::AbstractArray{<:Gray},
    algorithm::LopezAcosta=LopezAcosta(),
)
    landmask_dilated, landmask_ = create_landmask(
        landmask, algorithm.landmask_structuring_element
    )

    @info "Building cloudmask"
    cloudmask = create_cloudmask(falsecolor)

    # 2. Intermediate images
    @info "Finding ice labels"
    ice_labels = IceFloeTracker.find_ice_labels(falsecolor, landmask_dilated)

    @info "Sharpening truecolor image"
    # a. apply imsharpen to truecolor image using non-dilated landmask
    sharpened_truecolor = IceFloeTracker.imsharpen(truecolor, landmask_)
    # b. apply imsharpen to sharpened truecolor img using dilated landmask
    sharpened_gray_truecolor = IceFloeTracker.imsharpen_gray(
        sharpened_truecolor, landmask_dilated
    )

    @info "Normalizing truecolor image"
    normalized_image = IceFloeTracker.normalize_image(
        sharpened_truecolor, sharpened_gray_truecolor, landmask_dilated
    )

    # Discriminate ice/water
    @info "Discriminating ice/water"
    ice_water_discrim = IceFloeTracker.discriminate_ice_water(
        falsecolor, normalized_image, copy(landmask_dilated), cloudmask
    )

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    segA = IceFloeTracker.segmentation_A(
        IceFloeTracker.segmented_ice_cloudmasking(ice_water_discrim, cloudmask, ice_labels)
    )

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = IceFloeTracker.segmentation_B(sharpened_gray_truecolor, cloudmask, segA)

    # Process watershed in parallel using Folds
    @info "Building watersheds"
    watersheds_segB = Folds.map(
        IceFloeTracker.watershed_ice_floes, [segB.not_ice_bit, segB.ice_intersect]
    )
    # reuse the memory allocated for the first watershed
    watersheds_segB[1] .= IceFloeTracker.watershed_product(watersheds_segB...)

    # segmentation_F
    @info "Segmenting floes part 3/3"
    segF = IceFloeTracker.segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB[1],
        ice_labels,
        cloudmask,
        landmask_dilated,
    )

    labeled_floes = label_components(segF)

    # TODO: return ImageSegmentation.jl-style results
    return labeled_floes
end

@kwdef struct LopezAcostaTiling <: IceFloeTrackerAlgorithm
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

function segment(
    truecolor::{<:AbstractArray{<:Union{AbstractRGB,TransparentRGB}}},
    falsecolor::{<:AbstractArray{<:Union{AbstractRGB,TransparentRGB}}},
    landmask::AbstractArray{<:Gray},
    algorithm::LopezAcostaTiling=LopezAcostaTiling(),
)
    landmask_dilated, landmask_ = create_landmask(
        landmask, algorithm.landmask_structuring_element
    )

    @info "Remove alpha channel if it exists"
    rgb_truecolor_img = RGB.(truecolor)
    rgb_falsecolor_img = RGB.(falsecolor)

    @info "Get tile coordinates"
    tiles = IceFloeTracker.get_tiles(
        rgb_truecolor_img; rblocks=algorithm.tile_rblocks, cblocks=algorithm.tile_cblocks
    )
    @debug tiles

    @info "Set ice labels thresholds"
    ice_labels_thresholds = (
        prelim_threshold=algorithm.ice_labels_prelim_threshold,
        band_7_threshold=algorithm.ice_labels_band_7_threshold,
        band_2_threshold=algorithm.ice_labels_band_2_threshold,
        ratio_lower=algorithm.ice_labels_ratio_lower,
        ratio_upper=algorithm.ice_labels_ratio_upper,
        use_uint8=true,
    )
    @debug ice_labels_thresholds

    @info "Set adaptive histogram parameters"
    adapthisteq_params = (
        white_threshold=algorithm.adapthisteq_white_threshold,
        entropy_threshold=algorithm.adapthisteq_entropy_threshold,
        white_fraction_threshold=algorithm.adapthisteq_white_fraction_threshold,
    )
    @debug adapthisteq_params

    @info "Set gamma parameters"
    adjust_gamma_params = (
        gamma=algorithm.gamma,
        gamma_factor=algorithm.gamma_factor,
        gamma_threshold=algorithm.gamma_threshold,
    )
    @debug adjust_gamma_params

    @info "Set structuring elements"
    # This isn't tunable in the underlying code right now, 
    # so just use the defaults
    structuring_elements = IceFloeTracker.structuring_elements
    @debug structuring_elements

    @info "Set unsharp mask params"
    unsharp_mask_params = (
        radius=algorithm.unsharp_mask_radius,
        amount=algorithm.unsharp_mask_amount,
        factor=algorithm.unsharp_mask_factor,
    )
    @debug unsharp_mask_params

    @info "Set brighten factor"
    @debug brighten_factor

    @info "Set preliminary ice masks params"
    prelim_icemask_params = (
        radius=algorithm.prelim_icemask_radius,
        amount=algorithm.prelim_icemask_amount,
        factor=algorithm.prelim_icemask_factor,
    )
    @debug prelim_icemask_params

    @info "Set ice masks params"
    ice_masks_params = (
        band_7_threshold=algorithm.icemask_band_7_threshold,
        band_2_threshold=algorithm.icemask_band_2_threshold,
        band_1_threshold=algorithm.icemask_band_1_threshold,
        band_7_threshold_relaxed=algorithm.icemask_band_7_threshold_relaxed,
        band_1_threshold_relaxed=algorithm.icemask_band_1_threshold_relaxed,
        possible_ice_threshold=algorithm.icemask_possible_ice_threshold,
        k=algorithm.icemask_n_clusters, # number of clusters for kmeans segmentation
        factor=255, # normalization factor to convert images to uint8
    )
    @debug ice_masks_params

    @info "Segment floes"
    segmented_floes = IceFloeTracker.preprocess_tiling(
        n0f8.(rgb_falsecolor_img),
        n0f8.(rgb_truecolor_img),
        (; dilated=.!landmask_dilated, non_dilated=landmask_),
        tiles,
        ice_labels_thresholds,
        adapthisteq_params,
        adjust_gamma_params,
        structuring_elements,
        unsharp_mask_params,
        ice_masks_params,
        prelim_icemask_params,
        algorithm.brighten_factor,
    )

    labeled_floes = label_components(segmented_floes)

    return labeled_floes
end
