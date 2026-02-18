
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

@testitem "ContrastLimitedAdaptiveHistogramEqualization basic case" begin
    using TestImages
    using Images: adjust_histogram
    img = testimage("cameraman")
    img_eq = adjust_histogram(img, ContrastLimitedAdaptiveHistogramEqualization())
    @test size(img) == size(img_eq)
    @test 125075 ≈ sum(img_eq) rtol = 0.0001
end

@testitem "ContrastLimitedAdaptiveHistogramEqualization different image sizes and aspect ratios" begin
    # Test that the function can handle different sizes of input image without crashing
    using TestImages
    using Images: adjust_histogram, imresize
    
    function clahe_works_with_size(img, rsize, csize, rblocks, cblocks)
        img_resized = imresize(img, (rsize, csize))
        img_eq = adjust_histogram(img_resized, ContrastLimitedAdaptiveHistogramEqualization(;rblocks,cblocks))
        valid = size(img_resized) == size(img_eq)
        return valid
    end

    img = testimage("cameraman")
    # Happy Path – square images
    @test clahe_works_with_size(img, 1024, 1024, 8, 8)
    @test clahe_works_with_size(img, 512, 512, 8, 8)
    @test clahe_works_with_size(img, 256, 256, 8, 8)
    @test clahe_works_with_size(img, 128, 128, 8, 8)
    @test clahe_works_with_size(img, 64, 64, 8, 8)
    @test clahe_works_with_size(img, 32, 32, 8, 8)
    
    # Happy Path – non-square images
    @test clahe_works_with_size(img, 100, 78, 8, 8)
    @test clahe_works_with_size(img, 100, 78, 16, 16)
    @test clahe_works_with_size(img, 32, 78, 8, 8)
    @test clahe_works_with_size(img, 1024, 64, 32, 4)
    @test clahe_works_with_size(img, 512, 64, 32, 4)
    @test clahe_works_with_size(img, 2048, 1024, 100, 100)
    @test clahe_works_with_size(img, 33, 12, 3, 5)

    # Edge cases – lower limit of block size (3x3)
    @test clahe_works_with_size(img, 9, 9, 3, 3)
    @test clahe_works_with_size(img, 6, 6, 2, 2)
    @test clahe_works_with_size(img, 3, 3, 1, 1)

    # Broken cases
    @test clahe_works_with_size(img, 2, 2, 1, 1) broken=true # 2x2 blocks are too small
end
    


    # img = testimage("cameraman")
    # sizes = [1024, 512, 256, 128, 64, 32, 100, 78, 50, 20]
    # blocks = [1, 2, 4, 8, 16]
    # for rsize in sizes, csize in sizes, rblocks in blocks, cblocks in blocks
    #     f = ContrastLimitedAdaptiveHistogramEqualization(rblocks=rblocks, cblocks=cblocks)
    #     img_resized = imresize(img, (rsize, csize))
    #     img_eq = adjust_histogram(img_resized, f)
    #     @test size(img_resized) == size(img_eq)
    # end
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
    @test 27_311_946 ≈ sum(to_uint8(true_color_eq[:, :, 1])) rtol = 0.003

    # Use custom tile size
    side_length = size(true_color_eq, 1) ÷ 8
    tiles = get_tiles(true_color_image, side_length)
    true_color_eq = conditional_histeq(true_color_image, clouds_red, tiles)
    @test 30_255_658 ≈ sum(to_uint8(true_color_eq[:, :, 1])) rtol = 0.003
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
