abstract type IceFloeSegmentationAlgorithm end

struct LopezAcosta2019 <: IceFloeSegmentationAlgorithm
    landmask_structuring_element::AbstractMatrix{Bool}
end

function LopezAcosta2019(; landmask_structuring_element=make_landmask_se())
    return LopezAcosta2019(landmask_structuring_element)
end

function (p::LopezAcosta2019)(
    truecolor::T, falsecolor::T, landmask::U; return_intermediate_results::Bool=false
) where {T<:AbstractMatrix{<:AbstractRGB},U<:AbstractMatrix}

    # Move these conversions down through the function as each step gets support for 
    # the full range of image formats
    truecolor_image = float64.(truecolor)
    falsecolor_image = float64.(falsecolor)
    landmask_image = float64.(landmask)

    @info "building landmask"
    landmask_imgs = create_landmask(landmask_image, p.landmask_structuring_element)

    @info "Building cloudmask"
    cloudmask = create_cloudmask(falsecolor_image)

    # 2. Intermediate images
    @info "Finding ice labels"
    ice_labels = find_ice_labels(falsecolor_image, landmask_imgs.dilated)

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
        segmented_ice_cloudmasking(ice_water_discrim, cloudmask, ice_labels)
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
    @info "Segmenting floes part 3/3"
    segF = segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB_product,
        ice_labels,
        cloudmask,
        landmask_imgs.dilated,
    )

    @info "Labeling floes"
    labels_map = label_components(segF)

    # Return the original truecolor image, segmented
    segments = SegmentedImage(truecolor, labels_map)

    if return_intermediate_results
        intermediate_results = Dict(
            :landmask_dilated => landmask_imgs.dilated,
            :landmask_non_dilated => landmask_imgs.non_dilated,
            :cloudmask => cloudmask,
            :ice_labels => ice_labels,
            :sharpened_truecolor_image => sharpened_truecolor_image,
            :sharpened_gray_truecolor_image => sharpened_gray_truecolor_image,
            :normalized_image => normalized_image,
            :ice_water_discrim => ice_water_discrim,
            :segA => segA,
            :segB => segB,
            :watersheds_segB_product => watersheds_segB_product,
            :segF => segF,
            :labels_map => labels_map,
            :segments => segments,
        )
        return (segments, intermediate_results)
    else
        return segments
    end
end
