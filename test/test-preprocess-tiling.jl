
using IceFloeTracker:
    adapthisteq_params,
    adjust_gamma_params,
    brighten_factor,
    ice_labels_thresholds,
    ice_masks_params,
    prelim_icemask_params,
    preprocess_tiling,
    structuring_elements,
    unsharp_mask_params

    using IceFloeTracker:
    _get_masks,
    _process_image_tiles,
    _regularize,
    adjustgamma,
    apply_landmask,
    branchbridge,
    fillholes!,
    get_brighten_mask,
    get_final,
    get_holes,
    get_ice_masks,
    get_nlabel,
    get_segment_mask,
    get_tiles,
    histeq,
    imadjust,
    imbinarize,
    imbrighten,
    imcomplement,
    imgradientmag,
    impose_minima,
    imregionalmin,
    kmeans_segmentation,
    label_components,
    reconstruct,
    rgb2gray,
    se_disk2,
    se_disk4,
    to_uint8,
    unsharp_mask,
    watershed2
using Images

region = (1016:3045, 1486:3715)
data_dir = joinpath(@__DIR__, "test_inputs")
true_color_image = load(joinpath(data_dir, "NE_Greenland_truecolor.2020162.aqua.250m.tiff"));
ref_image = load(joinpath(data_dir, "NE_Greenland_reflectance.2020162.aqua.250m.tiff"));
landmask = float64.(load(joinpath(data_dir, "matlab_landmask.png"))) .> 0

true_color_image, ref_image, landmask = [
    img[region...] for img in (true_color_image, ref_image, landmask)
]

landmask = (dilated=landmask,)
tiles = get_tiles(true_color_image, rblocks=2, cblocks=3)

foo = preprocess_tiling(
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

@info "Step 1/2: Get masks"
begin
    mask_cloud_ice, clouds_view = _get_masks(float64.(ref_image); ice_labels_thresholds...)
    clouds_view .= .!mask_cloud_ice .* clouds_view

    # Get clouds_red for adaptive histogram equalization
    ref_img_cloudmasked = ref_image .* .!clouds_view
end

@info "Step 3: Tiled adaptive histogram equalization"
clouds_red = to_uint8(float64.(red.(ref_img_cloudmasked) .* 255))
clouds_red[.!landmask.dilated] .= 0
rgbchannels = _process_image_tiles(
    true_color_image, clouds_red, tiles, adapthisteq_params...
)
gammagreen = @view rgbchannels[:,:,2]
equalized_gray = rgb2gray(rgbchannels)


@info "Step 4: Remove clouds from equalized_gray"
masks = [f.(ref_img_cloudmasked) .== 0 for f in [red, green, blue]]
combo_mask = reduce((a, b) -> a .& b, masks)
equalized_gray[combo_mask] .= 0

@info "Step 5: unsharp_mask on equalized_gray and reconstruct"
sharpened = to_uint8(unsharp_mask(equalized_gray, unsharp_mask_params...))
equalized_gray_sharpened_reconstructed = reconstruct(
    sharpened, structuring_elements.se_disk1, "dilation", true
)
equalized_gray_sharpened_reconstructed[.!landmask.dilated] .= 0

# TODO: Steps 6 and 7 can be done in parallel as they are independent
@info "# Step 6: Repeat step 5 with equalized_gray"
equalized_gray_reconstructed = deepcopy(equalized_gray)
equalized_gray_reconstructed[.!landmask.dilated] .= 0
equalized_gray_reconstructed = reconstruct(
    equalized_gray_reconstructed, structuring_elements.se_disk1, "dilation", true
)
equalized_gray_reconstructed[.!landmask.dilated] .= 0

@info "STEP 7: Brighten equalized_gray"
brighten = get_brighten_mask(equalized_gray_reconstructed, gammagreen)
equalized_gray[.!landmask.dilated] .= 0
equalized_gray .= imbrighten(equalized_gray, brighten, brighten_factor)

@info "STEP 8: Get morphed_residue and adjust its gamma"
morphed_residue = clamp.(equalized_gray - equalized_gray_reconstructed, 0, 255)
agp = adjust_gamma_params
equalized_gray_sharpened_reconstructed_adjusted = imcomplement(
    adjustgamma(equalized_gray_sharpened_reconstructed, agp.gamma)
)
adjusting_mask = equalized_gray_sharpened_reconstructed_adjusted .> agp.gamma_threshold
morphed_residue[adjusting_mask] .=
    to_uint8.(morphed_residue[adjusting_mask] .* agp.gamma_factor)

@info "# Step 9: Get prelimnary ice masks"
prelim_icemask, binarized_tiling = get_ice_masks(
    ref_image, morphed_residue, landmask.dilated, tiles, true; ice_masks_params...
)

@info "Step 10: Get segmentation mask from preliminary icemask"
segment_mask = get_segment_mask(prelim_icemask, binarized_tiling)

@info "Step 11: Get local_maxima_mask and L0mask via watershed"
local_maxima_mask, L0mask = watershed2(morphed_residue, segment_mask, prelim_icemask)

@info "Step 12: Build icemask from all others"
local_maxima_mask = to_uint8(local_maxima_mask * 255)
prelim_icemask2 = get_combined_new(
    morphed_residue,
    local_maxima_mask,
    segment_mask,
    L0mask,
    structuring_elements.se_disk1;
    prelim_icemask_params...,
)

@info "Step 13: Get improved icemask"
icemask, _ = get_ice_masks(
    ref_image, prelim_icemask2, landmask.dilated, tiles, false; ice_masks_params...
)

@info "Step 14: Get final mask"
se = structuring_elements
se_erosion = se.se_disk1
se_dilation = se.se_disk2
final = get_final(icemask, segment_mask, se_erosion, se_dilation)

# return final