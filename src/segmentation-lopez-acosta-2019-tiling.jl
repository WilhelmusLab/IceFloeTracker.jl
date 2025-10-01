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

unsharp_mask_params = (radius=10, amount=2.0, factor=255.0)

brighten_factor = 0.1

ice_masks_params = (
    band_7_threshold=5 / 255,
    band_2_threshold=230 / 255,
    band_1_threshold=240 / 255,
    band_7_threshold_relaxed=10 / 255,
    band_1_threshold_relaxed=190 / 255,
    possible_ice_threshold=75 / 255,
    k=3, # number of clusters for kmeans segmentation
)

prelim_icemask_params = (radius=10, amount=2, factor=0.5)

diffusion_params = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")

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
    diffusion_params:AbstractDiffusionAlgorithm = diffusion_params
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

    begin
        @debug "Step 1/2: Create and apply cloudmask to reference image"

        cloudmask = IceFloeTracker.create_cloudmask(
            ref_image, LopezAcostaCloudMask(cloud_mask_thresholds...)
        )
        ref_img_cloudmasked = IceFloeTracker.apply_cloudmask(ref_image, cloudmask)
    end

    begin
        @debug "Step 3: Tiled adaptive histogram equalization"
        clouds_red = to_uint8(float64.(red.(ref_img_cloudmasked) .* 255))
        clouds_red[_landmask.dilated] .= 0

        # Apply Perona-Malik diffusion to each channel of true color image 
        # using the default inverse quadratic flux coefficient function
        true_color_diffused = IceFloeTracker.nonlinear_diffusion(
            float64.(true_color_image), p.diffusion_params
        )

        rgbchannels = conditional_histeq(
            true_color_diffused, clouds_red, tiles; adapthisteq_params...
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
        equalized_gray_sharpened_reconstructed[_landmask.dilated] .= 0
    end

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

    begin
        @debug "Step 8: Get morphed_residue and adjust its gamma"
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
        binarized_tiling = tiled_adaptive_binarization(Gray.(morphed_residue ./ 255), tiles;
                                           minimum_window_size=32, threshold_percentage=15) .> 0
        prelim_icemask = get_ice_masks(
            ref_image,
            Gray.(morphed_residue / 255),
            _landmask.dilated,
            tiles;
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
        icemask = get_ice_masks(
            ref_image,
            Gray.(prelim_icemask2 ./ 255),
            _landmask.dilated,
            tiles;
            ice_masks_params...,
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
