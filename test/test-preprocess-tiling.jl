

using IceFloeTracker: preprocess_tiling, ice_labels_thresholds, adapthisteq_params, adjust_gamma_params, structuring_elements, unsharp_mask_params, ice_masks_params, prelim_icemask_params, brighten_factor

region = (1016:3045, 1486:3715)
data_dir = joinpath(@__DIR__,"test_inputs")
tc_img = load(joinpath(data_dir,"NE_Greenland_truecolor.2020162.aqua.250m.tiff"))
ref_img = load(joinpath(data_dir,"NE_Greenland_reflectance.2020162.aqua.250m.tiff"))
landmask = load(joinpath(data_dir,"matlab_landmask.png"))

tc_img, ref_img, landmask = [img[region...] for img in (tc_img, ref_img, landmask)]

tc_img = tc_img[region...]

foo = preprocess_tiling(
    ref_img,
    tc_img,
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