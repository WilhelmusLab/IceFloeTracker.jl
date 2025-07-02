using Images
using IceFloeTracker:
    get_tiles,
    _process_image_tiles,
    to_uint8,
    unsharp_mask,
    imbrighten,
    get_ice_masks,
    imcomplement,
    adjustgamma,
    get_segment_mask,
    se_disk4,
    se_disk2,
    get_final,
    get_brighten_mask,
    reconstruct,
    histeq,
    label_components,
    watershed2,
    _regularize

@kwdef struct LopezAcosta2019Tiling <: IceFloeSegmentationAlgorithm
    # Landmask parameters
    landmask_structuring_element::AbstractMatrix{Bool} = make_landmask_se()

    structuring_elements = (
        se_disk1=collect(strel_diamond((3, 3))), se_disk2=se_disk2(), se_disk4=se_disk4()
    )

    # Tiling parameters
    tile_rblocks::Integer = 8
    tile_cblocks::Integer = 8

    # Ice labels thresholds
    ice_labels_prelim_threshold::Float64 = 110.0
    ice_labels_band_7_threshold::Float64 = 200.0
    ice_labels_band_2_threshold::Float64 = 190.0
    ice_labels_ratio_lower::Float64 = 0.0
    ice_labels_ratio_upper::Float64 = 0.75
    r_offset::Float64 = 0.0

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

function (p::LopezAcosta2019Tiling)(
    truecolor_image::T, falsecolor_image::T, landmask_image::U
) where {T<:Matrix{RGB{N0f8}},U<:AbstractMatrix}
    @info "Remove alpha channel if it exists"
    rgb_truecolor_img = RGB.(truecolor_image)
    rgb_falsecolor_img = RGB.(falsecolor_image)

    # Invert the landmasks – in the tiling version of the code, 
    # the landmask is expected to be the other polarity compared with
    # the non-tiling version.
    @info "building landmask"
    landmask_imgs = create_landmask(landmask_image, p.landmask_structuring_element)
    landmask = (dilated=.!landmask_imgs.dilated,)

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
        r_offset=p.r_offset,
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
    structuring_elements = p.structuring_elements
    @debug structuring_elements

    @info "Set unsharp mask params"
    unsharp_mask_params = (
        radius=p.unsharp_mask_radius,
        amount=p.unsharp_mask_amount,
        factor=p.unsharp_mask_factor,
    )
    @debug unsharp_mask_params

    @info "Set brighten factor"
    brighten_factor = p.brighten_factor
    @debug brighten_factor

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
        band_7_threshold_relaxed=p.icemask_band_7_threshold_relaxed,
        band_1_threshold_relaxed=p.icemask_band_1_threshold_relaxed,
        possible_ice_threshold=p.icemask_possible_ice_threshold,
        k=p.icemask_n_clusters, # number of clusters for kmeans segmentation
        factor=255, # normalization factor to convert images to uint8
    )
    @debug ice_masks_params

    @info "Segment floes"
    begin
        @debug "Step 1/2: Create and apply cloudmask to reference image"
        cloudmask = IceFloeTracker.create_cloudmask(
            rgb_falsecolor_img; ice_labels_thresholds...
        )
        ref_img_cloudmasked = IceFloeTracker.apply_cloudmask(rgb_falsecolor_img, cloudmask)
    end

    begin
        @debug "Step 3: Tiled adaptive histogram equalization"
        clouds_red = to_uint8(float64.(red.(ref_img_cloudmasked) .* 255))
        clouds_red[landmask.dilated] .= 0

        rgbchannels = _process_image_tiles(
            rgb_truecolor_img, clouds_red, tiles, adapthisteq_params...
        )

        gammagreen = @view rgbchannels[:, :, 2]
        equalized_gray = rgb2gray(rgbchannels)
    end

    begin
        @debug "Step 4: Remove clouds from equalized_gray"
        masks = [f.(ref_img_cloudmasked) .== 0 for f in [red, green, blue]]
        combo_mask = reduce((a, b) -> a .& b, masks)
        equalized_gray[combo_mask] .= 0
    end

    begin
        @debug "Step 5: unsharp_mask on equalized_gray and reconstruct"
        sharpened = to_uint8(unsharp_mask(equalized_gray, unsharp_mask_params...))
        equalized_gray_sharpened_reconstructed = reconstruct(
            sharpened, structuring_elements.se_disk1, "dilation", true
        )
        equalized_gray_sharpened_reconstructed[landmask.dilated] .= 0
    end

    # TODO: Steps 6 and 7 can be done in parallel as they are independent
    begin
        @debug "Step 6: Repeat step 5 with equalized_gray (landmasking, no sharpening)"
        equalized_gray_reconstructed = deepcopy(equalized_gray)
        equalized_gray_reconstructed[landmask.dilated] .= 0
        equalized_gray_reconstructed = reconstruct(
            equalized_gray_reconstructed, structuring_elements.se_disk1, "dilation", true
        )
        equalized_gray_reconstructed[landmask.dilated] .= 0
    end

    begin
        @debug "STEP 7: Brighten equalized_gray"
        brighten = get_brighten_mask(equalized_gray_reconstructed, gammagreen)
        equalized_gray[landmask.dilated] .= 0
        equalized_gray .= imbrighten(equalized_gray, brighten, brighten_factor)
    end

    begin
        @debug "STEP 8: Get morphed_residue and adjust its gamma"
        morphed_residue = clamp.(equalized_gray - equalized_gray_reconstructed, 0, 255)

        agp = adjust_gamma_params
        equalized_gray_sharpened_reconstructed_adjusted = imcomplement(
            adjustgamma(equalized_gray_sharpened_reconstructed, agp.gamma)
        )
        adjusting_mask =
            equalized_gray_sharpened_reconstructed_adjusted .> agp.gamma_threshold
        morphed_residue[adjusting_mask] .=
            to_uint8.(morphed_residue[adjusting_mask] .* agp.gamma_factor)
    end

    begin
        @debug "Step 9: Get preliminary ice masks"
        prelim_icemask, binarized_tiling = get_ice_masks(
            rgb_falsecolor_img,
            morphed_residue,
            landmask.dilated,
            tiles,
            true;
            ice_masks_params...,
        )
    end

    begin
        @debug "Step 10: Get segmentation mask from preliminary icemask"
        # Fill holes function in get_segment_mask a bit more aggressive than Matlabs
        segment_mask = get_segment_mask(prelim_icemask, binarized_tiling)
    end

    begin
        @debug "Step 11: Get local_maxima_mask and L0mask via watershed"
        local_maxima_mask, L0mask = watershed2(
            morphed_residue, segment_mask, prelim_icemask
        )
    end

    begin
        @debug "Step 12: Build icemask from all others"
        local_maxima_mask = to_uint8(local_maxima_mask * 255)
        prelim_icemask2 = _regularize(
            morphed_residue,
            local_maxima_mask,
            segment_mask,
            L0mask,
            structuring_elements.se_disk1;
            prelim_icemask_params...,
        )
    end

    begin
        @debug "Step 13: Get improved icemask"
        icemask, _ = get_ice_masks(
            rgb_falsecolor_img,
            prelim_icemask2,
            landmask.dilated,
            tiles,
            false;
            ice_masks_params...,
        )
    end

    begin
        @debug "Step 14: Get final mask"
        se = structuring_elements
        se_erosion = se.se_disk1
        se_dilation = se.se_disk2
        segmented_floes = get_final(icemask, segment_mask, se_erosion, se_dilation)
    end

    @info "Label floes"
    labeled_floes = label_components(segmented_floes)

    # TODO: return ImageSegmentation.jl-style results
    return (; labeled_floes, segmented_floes)
end
