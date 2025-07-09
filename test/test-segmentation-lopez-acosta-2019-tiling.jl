
using IceFloeTracker:
    adapthisteq_params,
    adjust_gamma_params,
    brighten_factor,
    cloud_mask_thresholds,
    ice_masks_params,
    prelim_icemask_params,
    preprocess_tiling,
    structuring_elements,
    unsharp_mask_params,
    get_tiles

using Images: labels_map, segment_mean

@testset "preprocess_tiling" begin
    region = (1016:3045, 1486:3714)
    data_dir = joinpath(@__DIR__, "test_inputs")
    true_color_image = load(
        joinpath(data_dir, "beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff")
    )
    ref_image = load(
        joinpath(data_dir, "beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff")
    )
    landmask = float64.(load(joinpath(data_dir, "matlab_landmask.png"))) .> 0

    # Crop images to region of interest
    true_color_image, ref_image, landmask = [
        img[region...] for img in (true_color_image, ref_image, landmask)
    ]

    landmask = (dilated=landmask,)
    tiles = get_tiles(true_color_image; rblocks=2, cblocks=3)

    segments = preprocess_tiling(
        ref_image,
        true_color_image,
        landmask,
        tiles,
        cloud_mask_thresholds,
        adapthisteq_params,
        adjust_gamma_params,
        structuring_elements,
        unsharp_mask_params,
        ice_masks_params,
        prelim_icemask_params,
        brighten_factor,
    )
    @info segments
    save(
        "./test_outputs/segmentation-Lopez-Acosta-2019-Tiling-detailed-segments-bitmatrix" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png",
        labels_map(segments) .> 0,
    )
    save(
        "./test_outputs/segmentation-Lopez-Acosta-2019-Tiling-detailed-segments-mean" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png",
        map(i -> segment_mean(segments, i), labels_map(segments)),
    )
    # dmw: replace with test of mismatch against a preprocessed image
    @test abs(sum(labels_map(segments) .> 0) - 1461116) / 1461116 < 0.1
end
