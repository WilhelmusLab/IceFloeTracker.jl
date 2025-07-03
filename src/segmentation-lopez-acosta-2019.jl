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
        landmask_imgs.dilated,
    )

    labels_map = label_components(segF)
    segments = ImageSegmentation.SegmentedImage(truecolor_image, labels_map)

    return segments
end