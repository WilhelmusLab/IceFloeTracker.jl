
@testitem "Create Cloudmask" begin
    using IceFloeTracker: @test_approx_eq_sigma_eps
    using Images: RGBA, N0f8
    include("config.jl")

    # define constants, maybe move to test config file
    matlab_cloudmask_file = "$(test_data_dir)/matlab_cloudmask.tiff"

    # Create and apply cloudmask
    ref_image = float64.(load(falsecolor_test_image_file)[test_region...])

    matlab_cloudmask = float64.(load(matlab_cloudmask_file))
    @time cloudmask = IceFloeTracker.create_cloudmask(ref_image)
    @time masked_image = IceFloeTracker.apply_cloudmask(
        ref_image, cloudmask; modify_channel_1=true
    )

    # test for percent difference in cloudmask images
    @test (@test_approx_eq_sigma_eps masked_image matlab_cloudmask [0, 0] 0.005) === nothing

    # test for create_clouds_channel
    clouds_channel_expected = load(clouds_channel_test_file)
    clds_channel = IceFloeTracker.create_clouds_channel(cloudmask, ref_image)
    @test (@test_approx_eq_sigma_eps (clds_channel) (clouds_channel_expected) [0, 0] 0.005) ===
        nothing

    @info "Test image that loads as RGBA"
    pth_RGBA_tiff = "$(test_data_dir)/466-sea_of_okhostk-100km-20040421.terra.truecolor.250m.tiff"
    ref_image = load(pth_RGBA_tiff)
    @test typeof(ref_image) <: Matrix{RGBA{N0f8}}
    cloudmask = IceFloeTracker.create_cloudmask(ref_image)
    @test sum(.!cloudmask) === 0 # all pixels are clouds
end

# Test creation and application of multiple cloudmask types
@testitem "Cloudmask Customization" begin
    using IceFloeTracker:
        Watkins2025GitHub, LopezAcostaCloudMask, Watkins2025CloudMask, create_cloudmask

    data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    case = first(data_loader(c -> (c.case_number == 6 && c.satellite == "terra")))

    # Settings from Watkins et al. 2025
    cloud_mask_settings = (
        prelim_threshold=53.0 / 255.0,
        band_7_threshold=130.0 / 255.0,
        band_2_threshold=170.0 / 255.0,
        ratio_lower=0.0,
        ratio_offset=0.0,
        ratio_upper=0.52,
    )

    cmask_orig = LopezAcostaCloudMask(cloud_mask_settings...)
    cmask_morpho = Watkins2025CloudMask()

    cmask_orig_img = create_cloudmask(case.modis_falsecolor, cmask_orig)
    cmask_morpho_img = create_cloudmask(case.modis_falsecolor, cmask_morpho)
    @test sum(cmask_orig_img) >= sum(cmask_morpho_img)
end
