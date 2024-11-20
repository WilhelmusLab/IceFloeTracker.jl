
using IceFloeTracker:
    adapthisteq_params,
    adjust_gamma_params,
    brighten_factor,
    ice_labels_thresholds,
    ice_masks_params,
    prelim_icemask_params,
    preprocess_tiling,
    structuring_elements,
    unsharp_mask_params,
    get_tiles

@testset "preprocess_tiling" begin
    region = (1016:3045, 1486:3714)
    data_dir = joinpath(@__DIR__, "test_inputs")
    true_color_image = load(
        joinpath(data_dir, "NE_Greenland_truecolor.2020162.aqua.250m.tiff")
    )
    ref_image = load(joinpath(data_dir, "NE_Greenland_reflectance.2020162.aqua.250m.tiff"))
    landmask = float64.(load(joinpath(data_dir, "matlab_landmask.png"))) .> 0

    # Crop images to region of interest
    true_color_image, ref_image, landmask = [
        img[region...] for img in (true_color_image, ref_image, landmask)
    ]

    landmask = (dilated=landmask,)
    tiles = get_tiles(true_color_image; rblocks=2, cblocks=3)

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

    @test sum(foo) == 1735472
end
