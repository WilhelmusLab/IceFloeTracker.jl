
@testsnippet FalseColorCloudmask begin
    import IceFloeTracker.Preprocessing: _get_masks
    using Images: channelview
    """
    Private function for testing the conditional adaptive histogram equalization workflow.
    """
    function _get_false_color_cloudmasked(;
        false_color_image,
        prelim_threshold=110.0,
        band_7_threshold=200.0,
        band_2_threshold=190.0,
    )
        mask_cloud_ice, clouds_view = _get_masks(
            false_color_image;
            prelim_threshold=prelim_threshold / 255.0,
            band_7_threshold=band_7_threshold / 255.0,
            band_2_threshold=band_2_threshold / 255.0,
            ratio_lower=0.0,
            ratio_offset=0.0,
            ratio_upper=0.75,
        )

        clouds_view[mask_cloud_ice] .= 0

        # remove clouds and land from each channel
        channels = Int.(channelview(false_color_image) * 255)

        # Apply the mask to each channel
        for i in 1:3
            @views channels[i, :, :][clouds_view] .= 0
        end

        return channels
    end
end

@testitem "conditional adaptivehisteq (data loader)" setup = [FalseColorCloudmask] begin
    dataset = filter(
        c -> c.case_number == 161 && c.satellite == "terra", Watkins2026Dataset()
    )
    case = first(dataset)
    false_color_image = modis_falsecolor(case)
    true_color_image = modis_truecolor(case)
    landmask = modis_landmask(case)

    clouds = _get_false_color_cloudmasked(;
        false_color_image=false_color_image,
        prelim_threshold=110.0,
        band_7_threshold=200.0,
        band_2_threshold=190.0,
    )
    clouds_red = clouds[1, :, :]

    dilated_landmask = create_landmask(landmask).dilated

    clouds_red[dilated_landmask] .= 0
    @test sum(clouds_red) ≈ 10_350_341 rtol = 0.01

    tiles = get_tiles(true_color_image; rblocks=2, cblocks=2)

    true_color_eq = conditional_histeq(true_color_image, clouds_red, tiles)

    # This differs from MATLAB script due to disparity in the implementations
    # of the adaptive histogram equalization / diffusion functions
    # For the moment testing for regression
    @test sum(to_uint8(true_color_eq[:, :, 1])) ≈ 27_422_448 rtol = 0.003

    # Use custom tile size
    side_length = size(true_color_eq, 1) ÷ 8
    tiles = get_tiles(true_color_image, side_length)
    true_color_eq = conditional_histeq(true_color_image, clouds_red, tiles)
    @test sum(to_uint8(true_color_eq[:, :, 1])) ≈ 27_446_614 rtol = 0.003
end

@testitem "_get_false_color_cloudmasked (data loader)" setup = [FalseColorCloudmask] begin
    dataset = filter(
        c -> c.case_number == 161 && c.satellite == "terra", Watkins2026Dataset()
    )
    case = first(dataset)
    false_color_image = modis_falsecolor(case)

    false_color_cloudmasked = _get_false_color_cloudmasked(;
        false_color_image=false_color_image,
        prelim_threshold=110.0,
        band_7_threshold=200.0,
        band_2_threshold=190.0,
    )

    @test sum(false_color_cloudmasked[1, :, :]) ≈ 10_350_341 rtol = 0.01
    @test sum(false_color_cloudmasked[2, :, :]) ≈ 17_029_014 rtol = 0.01
    @test sum(false_color_cloudmasked[3, :, :]) ≈ 17_159_247 rtol = 0.01
end

@testitem "adapthisteq" begin
    using TestImages: testimage
    function convert_to_255_matrix(img)::Matrix{Int}
        img_clamped = clamp.(img, 0.0, 1.0)
        return round.(Int, img_clamped * 255)
    end

    img = convert_to_255_matrix(testimage("cameraman"))
    img_eq = adapthisteq(img)
    @test sum(img_eq) ≈ 32_387_397 rtol = 0.05
end

@testitem "rgb2gray" begin
    using Images: RGB, N0f8, float64
    pixels = float64.([RGB{N0f8}(0.345, 0.361, 0.404) RGB{N0f8}(0.808, 0.808, 0.8);])
    h = rgb2gray(pixels)
    @test h[1, 1] == 92
    @test h[1, 2] == 206
end

@testitem "histeq" begin
    @testset "normal cases" begin
        @test histeq([
            4 4 4 4 4
            3 4 5 4 3
            3 5 5 5 3
            3 4 5 4 3
            4 4 4 4 4
        ]) == [
            204 204 204 204 204
            061 204 255 204 061
            061 255 255 255 061
            061 204 255 204 061
            204 204 204 204 204
        ]
    end

    @testset "edge cases" begin
        @test histeq([
            0 0 0
            2 2 2
        ]) == [
            128 128 128
            255 255 255
        ]

        @test histeq([
            000 000 000
            255 255 255
        ]) == [
            128 128 128
            255 255 255
        ]
    end
end
