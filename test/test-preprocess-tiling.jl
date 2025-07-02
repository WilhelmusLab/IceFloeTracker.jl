
using IceFloeTracker: LopezAcosta2019Tiling

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

    algorithm = LopezAcosta2019Tiling(;
        tile_rblocks=2,
        tile_cblocks=3,

        # Ice labels thresholds
        ice_labels_prelim_threshold=110.0,
        ice_labels_band_7_threshold=200.0,
        ice_labels_band_2_threshold=190.0,
        ice_labels_ratio_lower=0.0,
        ice_labels_ratio_upper=0.75,
        r_offset=0.0,

        # Adaptive histogram equalization parameters
        adapthisteq_white_threshold=25.5,
        adapthisteq_entropy_threshold=4,
        adapthisteq_white_fraction_threshold=0.4,

        # Gamma parameters,
        gamma=1.5,
        gamma_factor=1.3,
        gamma_threshold=220.0,

        # Unsharp mask parameters,
        unsharp_mask_radius=10,
        unsharp_mask_amount=2.0,
        unsharp_mask_factor=255.0,

        # Brighten parameters,
        brighten_factor=0.1,

        # Preliminary ice mask parameters,
        prelim_icemask_radius=10,
        prelim_icemask_amount=2,
        prelim_icemask_factor=0.5,

        # Main ice mask parameters,
        icemask_band_7_threshold=5,
        icemask_band_2_threshold=230,
        icemask_band_1_threshold=240,
        icemask_band_7_threshold_relaxed=10,
        icemask_band_1_threshold_relaxed=190,
        icemask_possible_ice_threshold=75,
        icemask_n_clusters=3,
    )
    (labeled_floes, segmented_floes) = algorithm(ref_image, true_color_image, landmask)

    @test sum(segmented_floes) == 1461116
end
