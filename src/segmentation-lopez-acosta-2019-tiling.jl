using Images
using IceFloeTracker:
    get_tiles,
    conditional_histeq,
    to_uint8,
    unsharp_mask,
    imbrighten,
    get_ice_masks,
    imcomplement,
    adjustgamma,
    to_uint8,
    get_segment_mask,
    se_disk4,
    se_disk2,
    branchbridge,
    fillholes!,
    get_final,
    apply_landmask,
    apply_cloudmask,
    kmeans_segmentation,
    get_brighten_mask,
    reconstruct,
    imgradientmag,
    histeq,
    label_components,
    watershed2,
    tiled_adaptive_binarization,
    _regularize

# Sample input parameters expected by the main function
cloud_mask_thresholds = (
    prelim_threshold=110.0 / 255.0,
    band_7_threshold=200.0 / 255.0,
    band_2_threshold=190.0 / 255.0,
    ratio_lower=0.0,
    ratio_offset=0.0,
    ratio_upper=0.75,
)

adapthisteq_params = (
    white_threshold=25.5, entropy_threshold=4, white_fraction_threshold=0.4
)

adjust_gamma_params = (gamma=1.5, gamma_factor=1.3, gamma_threshold=220)

structuring_elements = (
    se_disk1=collect(strel_diamond((3, 3))), se_disk2=se_disk2(), se_disk4=se_disk4()
)

unsharp_mask_params = (radius=10, amount=2.0, threshold=0.0)

brighten_factor = 0.1 # brightens by darkening the complement

ice_masks_params = (
    band_7_threshold=5/255,
    band_2_threshold=230/255,
    band_1_threshold=240/255,
    band_7_threshold_relaxed=10/255,
    band_1_threshold_relaxed=190/255,
    possible_ice_threshold=75/255
)

prelim_icemask_params = (radius=10, amount=2, factor=[0.3, 0.5])

@kwdef struct LopezAcosta2019Tiling <: IceFloeSegmentationAlgorithm
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

function (p::LopezAcosta2019Tiling)(
    truecolor::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    falsecolor::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    landmask::AbstractArray{<:Union{AbstractGray,AbstractRGB,TransparentRGB}};
    intermediate_results_callback::Union{Nothing,Function}=nothing,
)
    _landmask = IceFloeTracker.create_landmask(landmask, strel_box((3,3))) # smaller strel than in some test cases
    tiles = get_tiles(truecolor; p.tile_settings...)

    ref_image = RGB.(falsecolor)  # TODO: remove this typecast
    true_color_image = RGB.(truecolor)  # TODO: remove this typecast
    true_color_diffused = IceFloeTracker.nonlinear_diffusion(float64.(true_color_image), 0.1, 75, 3)

    begin
        @debug "Step 1/2: Create and apply cloudmask to reference image"

        cloudmask = IceFloeTracker.create_cloudmask(
            ref_image, LopezAcostaCloudMask(cloud_mask_thresholds...)
        )
        ref_img_cloudmasked = apply_landmask(
                                apply_cloudmask(ref_image, cloudmask),
                                .!_landmask.dilated) # TODO: Clarify landmask 1/0 meaning
        
    end

    begin
        @debug "Step 3: Tiled adaptive histogram equalization"
        rgbchannels = channelview(true_color_diffused)
        # tiles = filter(test_function, tiles) Add a test_function that selects only the tiles with ocean pixels
        
        for tile in tiles
            clouds_tile = red.(ref_img_cloudmasked[tile...])
            entropy = Images.entropy(clouds_tile) # Entropy calculation works on grayscale images and on uint8
            whitefraction = sum(clouds_tile .> adapthisteq_params.white_threshold / 255) / length(clouds_tile) # threshold depends on image type. Update to 0-1.

        rgbchannels = conditional_histeq(
            true_color_image, clouds_red, tiles; adapthisteq_params...
        )

        true_color_equalized = colorview(eltype(true_color_diffused), rgbchannels)

        gammagreen = green.(true_color_equalized)
        equalized_gray = Gray.(true_color_equalized)
    end

    # dmw: Replacing map reduction step. Note, I don't think this step is necessary here. Maybe at the end.
    begin
        @debug "Step 4: Apply cloudmask and landmask to the equalized image"
        apply_cloudmask!(equalized_gray, cloudmask)
        apply_landmask!(equalized_gray, .!_landmask.dilated)
    end
   
    # begin
    #     @debug "Step 4: Remove clouds from equalized_gray"
    #     masks = [f.(ref_img_cloudmasked) .== 0 for f in [red, green, blue]]
    #     combo_mask = reduce((a, b) -> a .& b, masks)
    #     equalized_gray[combo_mask] .= 0
    # end

    # TODO: Steps 5 and 6 can be done in parallel as they are independent
    begin
        @debug "Step 5: Reconstruct equalized gray by dilation"
        dilated_img = dilate(equalized_gray, strel_diamond((3, 3)))
        equalized_gray_reconstructed = mreconstruct(dilate, complement.(dilated_img), complement.(equalized_gray), strel_diamond((3, 3)))
        apply_landmask!(equalized_gray_reconstructed, .!_landmask.dilated)
    end

    begin
        @debug "Step 6: unsharp_mask on equalized_gray and reconstruct by dilation"
        sharpened_img = unsharp_mask(equalized_gray, unsharp_mask_params...)
        dilated_img = dilate(sharpened_img, strel_diamond((3, 3)))
        equalized_gray_sharpened_reconstructed = mreconstruct(dilate, complement.(dilated_img), complement.(sharpened_img), strel_diamond((3, 3)))
        apply_landmask!(equalized_gray_sharpened_reconstructed, .!_landmask.dilated)
    end

    begin
        @debug "Step 7: Brighten equalized_gray and compute residue"
        # dmw: these steps are too simple to have dedicated functions imo
        # further, it's confusing that we are brightening by darkening a complement. 
        # the mask selects where the pixels in the reconstruction (dark floes, bright leads)
        # are brighter than the original image. so what is happening is that leads and gaps between
        # floes are being darkened by a multiplicative factor "bright factor".
        # Note: very sensitive to the data being float, not N0f8

    # TODO: Steps 6 and 7 can be done in parallel as they are independent
    begin
        @debug "Step 6: Repeat step 5 with equalized_gray (landmasking, no sharpening)"
        equalized_gray_reconstructed = deepcopy(equalized_gray)
        IceFloeTracker.apply_landmask!(equalized_gray_reconstructed, _landmask.dilated)

        equalized_gray_reconstructed = reconstruct(
            equalized_gray_reconstructed, structuring_elements.se_disk1, "dilation", true
        )
        IceFloeTracker.apply_landmask!(equalized_gray_reconstructed, _landmask.dilated)
    end

    begin
        @debug "Step 7: Brighten equalized_gray"
        brighten = get_brighten_mask(equalized_gray_reconstructed, gammagreen)
        IceFloeTracker.apply_landmask!(equalized_gray, _landmask.dilated)
        equalized_gray .= imbrighten(equalized_gray, brighten, brighten_factor)
    end

        equalized_gray_sharpened_reconstructed_adjusted = complement.(
            adjust_histogram(equalized_gray_sharpened_reconstructed, GammaCorrection(agp.gamma)))
        adjusting_mask =
            equalized_gray_sharpened_reconstructed_adjusted .> agp.gamma_threshold ./ 255 # gamma threshold depends on image type
        morphed_residue[adjusting_mask] .= morphed_residue[adjusting_mask] .* agp.gamma_factor
        clamp!(morphed_residue, 0, 1)
        # morphed_residue .= Gray.(morphed_residue ./ maximum(morphed_residue))
    end

    begin
        @debug "Step 9: Get preliminary ice masks"
        binarized_tiling = tiled_adaptive_binarization(Gray.(morphed_residue ./ 255), tiles;
                                           minimum_window_size=32, threshold_percentage=15) .> 0
        prelim_icemask = get_ice_masks(
            ref_image, Gray.(morphed_residue), _landmask.dilated, tiles; ice_masks_params..., k=4
        )
    end

    begin
        @debug "Step 10: Get segmentation mask from preliminary icemask"
        # Fill holes function in get_segment_mask a bit more aggressive 
        segment_mask = get_segment_mask(prelim_icemask, binarized_tiling)
    end

    begin # _reconst_watershed requires an integer matrix
        @debug "Step 11: Get local_maxima_mask and L0mask via watershed"
        local_maxima_mask, L0mask = watershed2(
            morphed_residue, segment_mask, prelim_icemask
        )
    end

    begin
        @debug "Step 12: Build icemask from all others"
        # local_maxima_mask = to_uint8(local_maxima_mask * 255) # dmw: lmm is binary, I think
        prelim_icemask2 = _regularize(
            morphed_residue, 
            local_maxima_mask,
            segment_mask,
            L0mask,
            structuring_elements.se_disk1;
            prelim_icemask_params..., # choice of factor? It's 0.3 in the MATLAB code.
        )
    end

    begin
        @debug "Step 13: Get improved icemask"
        icemask = get_ice_masks(
            ref_image, Gray.(prelim_icemask2), _landmask.dilated, tiles; ice_masks_params..., k=3
        )
    end

    begin
        @debug "Step 14: Get final mask"
        se = structuring_elements
        se_erosion = se.se_disk1
        se_dilation = se.se_disk2
        binary_floe_masks = get_final(icemask, segment_mask, se_erosion, se_dilation)
    end

    begin
        @debug "Step 15: Converting mask to segmented image"
        labels = label_components(binary_floe_masks)
        segmented = SegmentedImage(truecolor, labels)
    end

    if !isnothing(intermediate_results_callback)
        segments_truecolor = SegmentedImage(truecolor, labels)
        segments_falsecolor = SegmentedImage(falsecolor, labels)
        intermediate_results_callback(;
            falsecolor,
            truecolor,
            ref_img_cloudmasked,
            gammagreen,
            equalized_gray,
            equalized_gray_sharpened_reconstructed,
            equalized_gray_reconstructed,
            morphed_residue,
            prelim_icemask,
            binarized_tiling,
            segment_mask,
            local_maxima_mask,
            L0mask,
            prelim_icemask2,
            icemask,
            binary_floe_masks,
            labels,
            segmented,
            segment_mean_truecolor=map(i -> segment_mean(segments_truecolor, i), labels),
            segment_mean_falsecolor=map(i -> segment_mean(segments_falsecolor, i), labels),
        )
    end

    return segmented
end
