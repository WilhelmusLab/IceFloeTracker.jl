
begin
    using Serialization
    using Images
    using IceFloeTracker:
        get_tiles,
        _get_masks,
        _process_image_tiles,
        to_uint8,
        unsharp_mask,
        reconstruct_erosion,
        imbrighten,
        imadjust,
        get_ice_masks,
        imcomplement,
        adjustgamma,
        to_uint8,
        get_holes,
        get_segment_mask,
        se_disk4

    tc_img = load("test/test_inputs/NE_Greenland_truecolor.2020162.aqua.250m.tiff")
    ref_img = deserialize("ref_image.jls")
    landmask = deserialize("landmasks.jls")
    tiles = IceFloeTracker.get_tiles(ref_img; rblocks=8, cblocks=6)

    # ice_labels_thresholds = (
    prelim_threshold = 110.0
    band_7_threshold = 200.0
    band_2_threshold = 190.0
    ratio_lower = 0.0
    ratio_upper = 0.75
    use_uint8 = true
    # )

    # adapthisteq params
    white_threshold = 25.5
    entropy_threshold = 4
    white_fraction_threshold = 0.4

    # adjustgamma params
    gamma = 1.5
    gamma_factor = 1.3
    gamma_threshold = 220
end

#= Preprocesing with tiling

Inputs

- ref_img: Reference image
- tc_img
- landmask: Landmask
- tiling
- ice_labels_thresholds
function preprocess_tiling(ref_image, true_color_image, landmask, ice_labels_thresholds)
=#
#= Get these two masks
clouds_view is masked with mask_cloud_ice =#
# Step 1/2
begin
    mask_cloud_ice, clouds_view = _get_masks(
        ref_img;
        prelim_threshold=prelim_threshold,
        band_7_threshold=band_7_threshold,
        band_2_threshold=band_2_threshold,
        ratio_lower=ratio_lower,
        ratio_upper=ratio_upper,
        use_uint8=use_uint8,
    )
    clouds_view .= .!mask_cloud_ice .* clouds_view
    @assert sum(mask_cloud_ice) == 23385158
    @assert 4177081 == sum(.!mask_cloud_ice .* clouds_view)

    # Step 3: Get clouds_red for adaptive histogram equalization
    ref_img_cloudmasked = ref_img .* .!clouds_view
    # channelview(ref_img_cloudmasked)[3, :, :] * 255 |> sum
end

# Step 4: Tiled adaptive histogram equalization
clouds_red = to_uint8(float64.(red.(ref_img_cloudmasked) .* 255))
clouds_red[.!landmask.dilated] .= 0

@time equalized_gray, gammagreen = _process_image_tiles(
    tc_img, clouds_red, tiles, white_threshold, entropy_threshold, white_fraction_threshold
)

# Step 5: Remove clouds from equalized_gray
begin
    masks = [f.(ref_img_cloudmasked) .== 0 for f in [red, green, blue]]
    combo_mask = reduce((a, b) -> a .& b, masks)
    @assert sum(combo_mask) == 5708073
    equalized_gray[combo_mask] .= 0
    sum(equalized_gray)
end

# Step 6: unsharp_mask on equalized_gray
begin
    equalized_gray_sharpened = to_uint8(
        IceFloeTracker.unsharp_mask(equalized_gray, 10, 2.0, 255)
    )
    Gray.(equalized_gray_sharpened / 255)
end

# Step 7: Apply reconstruct_erosion on equalized_gray_sharpened and landmask (normalization)
begin
    se_disk1 = IceFloeTracker.MorphSE.StructuringElements.strel_diamond((3, 3))
    se = se_disk1
    equalized_gray_sharpened_reconstructed = IceFloeTracker.reconstruct_erosion(
        equalized_gray_sharpened, se_disk1
    )
    equalized_gray_sharpened_reconstructed[.!landmask.dilated] .= 0
end

# Step 8: Repeat step 7 with equalized_gray
begin
    equalized_gray_reconstructed = deepcopy(equalized_gray)
    equalized_gray_reconstructed[.!landmask.dilated] .= 0
    @time equalized_gray_reconstructed = IceFloeTracker.reconstruct_erosion(
        equalized_gray_reconstructed, se_disk1
    )
    equalized_gray_reconstructed[.!landmask.dilated] .= 0
end

# STEP 9: Get brighten
brighten = equalized_gray_reconstructed - gammagreen

# STEP 10: Get equalized_gray_bright
begin
    equalized_gray[.!landmask.dilated] .= 0
    equalized_gray .= imbrighten(equalized_gray, brighten, 0.1)
end

# STEP 11: Get morphed_residue
morphed_residue = clamp.(equalized_gray - equalized_gray_reconstructed, 0, 255)

# STEP 12: Adjust gamma on morphed_residue
begin
    equalized_gray_sharpened_reconstructed_adjusted = imcomplement(
        adjustgamma(equalized_gray_sharpened_reconstructed, gamma)
    )
    adjusting_mask = equalized_gray_sharpened_reconstructed .> gamma_threshold
    morphed_residue[adjusting_mask] .=
        to_uint8.(morphed_residue[adjusting_mask] .* gamma_factor)
end

# Step 13: Get ice masks
begin
    binarize::Bool = true
    band_7_threshold = 5
    band_2_threshold = 230
    band_1_threshold = 240
    band_7_threshold_relaxed = 10
    band_1_threshold_relaxed = 190
    possible_ice_threshold = 75
    factor = 255

    @time prelim_icemask, binarized_tiling, _ = get_ice_masks(
        ref_img,
        morphed_residue,
        landmask.dilated,
        tiles,
        binarize;
        band_7_threshold=band_7_threshold,
        band_2_threshold=band_2_threshold,
        band_1_threshold=band_1_threshold,
        band_7_threshold_relaxed=band_7_threshold_relaxed,
        band_1_threshold_relaxed=band_1_threshold_relaxed,
        possible_ice_threshold=possible_ice_threshold,
        factor=factor,
    )
end

# Step 14: Get segment mask
segment_mask = get_segment_mask(prelim_icemask, binarized_tiling)

# Step 15: watershed
fgm, L0mask = IceFloeTracker.watershed2(morphed_residue, segment_mask)

# Step 16: Get new3
radius = 10
amount = 2
factor = .5

fgm = to_uint8(fgm * 255)
new2 = get_new2(morph_residue, fgm, factor, segment_mask, L0mask)
new3 = get_new3(new2, L0mask, radius, amount, local_maxima_mask, factor, segment_mask, se_disk1)

# Step 17: Get improved icemask
@time icemask, _, label = get_ice_masks(
    ref_img,
    new3,
    landmask.dilated,
    tiles,
    false;
    band_7_threshold=band_7_threshold,
    band_2_threshold=band_2_threshold,
    band_1_threshold=band_1_threshold,
    band_7_threshold_relaxed=band_7_threshold_relaxed,
    band_1_threshold_relaxed=band_1_threshold_relaxed,
    possible_ice_threshold=possible_ice_threshold,
    factor=factor,
)

# Step 18: Get final
se_erosion = se_disk1
se_dilation = se_disk4()
final = get_final(icemask, label, segment_mask, se_erosion, se_dilation)