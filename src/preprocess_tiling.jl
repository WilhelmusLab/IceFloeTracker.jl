using Images
using IceFloeTracker:
    get_tiles,
    _get_masks,
    _process_image_tiles,
    to_uint8,
    unsharp_mask,
    imbrighten,
    imadjust,
    get_ice_masks,
    imcomplement,
    adjustgamma,
    to_uint8,
    get_holes,
    get_segment_mask,
    se_disk4,
    se_disk2,
    branchbridge,
    fillholes!,
    get_final,
    apply_landmask,
    kmeans_segmentation,
    get_nlabel,
    get_brighten_mask,
    get_holes,
    reconstruct,
    imgradientmag,
    histeq,
    impose_minima,
    label_components,
    imregionalmin,
    watershed2,
    imbinarize,
    _regularize

# Sample input parameters expected by the main function
ice_labels_thresholds = (
    prelim_threshold=110.0,
    band_7_threshold=200.0,
    band_2_threshold=190.0,
    ratio_lower=0.0,
    ratio_upper=0.75,
    r_offset=0.0,
)

adapthisteq_params = (
    white_threshold=25.5, entropy_threshold=4, white_fraction_threshold=0.4
)

adjust_gamma_params = (gamma=1.5, gamma_factor=1.3, gamma_threshold=220)

structuring_elements = (
    se_disk1=collect(strel_diamond((3, 3))), se_disk2=se_disk2(), se_disk4=se_disk4()
)

unsharp_mask_params = (radius=10, amount=2.0, factor=255.0)

brighten_factor = 0.1

ice_masks_params = (
    band_7_threshold=5,
    band_2_threshold=230,
    band_1_threshold=240,
    band_7_threshold_relaxed=10,
    band_1_threshold_relaxed=190,
    possible_ice_threshold=75,
    k=3, # number of clusters for kmeans segmentation
    factor=255, # normalization factor to convert images to uint8
)

prelim_icemask_params = (radius=10, amount=2, factor=0.5)

function preprocess_tiling(
    ref_image,
    true_color_image,
    landmask,
    tiles,
    ice_labels_thresholds,
    adapthisteq_params,
    adjust_gamma_params,
    structuring_elements,
    unsharp_mask_params,
    ice_masks_params,
    prelim_icemask_params,
    brighten_factor,
)
    begin
        @debug "Step 1/2: Create and apply cloudmask to reference image"
        cloudmask = IceFloeTracker.create_cloudmask(ref_image; ice_labels_thresholds...)
        ref_img_cloudmasked = IceFloeTracker.apply_cloudmask(ref_image, cloudmask)
    end

    begin
        @debug "Step 3: Tiled adaptive histogram equalization"
        clouds_red = to_uint8(float64.(red.(ref_img_cloudmasked) .* 255))
        clouds_red[landmask.dilated] .= 0

        rgbchannels = _process_image_tiles(
            true_color_image, clouds_red, tiles, adapthisteq_params...
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
            ref_image, morphed_residue, landmask.dilated, tiles, true; ice_masks_params...
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
            ref_image, prelim_icemask2, landmask.dilated, tiles, false; ice_masks_params...
        )
    end

    begin
        @debug "Step 14: Get final mask"
        se = structuring_elements
        se_erosion = se.se_disk1
        se_dilation = se.se_disk2
        final = get_final(icemask, segment_mask, se_erosion, se_dilation)
    end

    return final
end
