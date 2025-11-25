
@testitem "Conditional adaptivehisteq" begin
    import IceFloeTracker.Preprocessing: _get_masks
    using Images: load, float64, channelview
    using TestImages: testimage

    include("config.jl")

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

    begin
        datadir = joinpath(@__DIR__, "test_inputs/")
        path_true_color_image = joinpath(
            datadir, "beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff"
        )
        path_false_color_image = joinpath(
            datadir, "beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff"
        )
        true_color_image = float64.(load(path_true_color_image))
        false_color_image = float64.(load(path_false_color_image))
        dilated_landmask = BitMatrix(load(joinpath(datadir, "matlab_landmask_dilated.png")))
    end

    @testset "Prereq cloud image" begin
        false_color_cloudmasked = _get_false_color_cloudmasked(;
            false_color_image=false_color_image,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0,
        )

        # replaced exact equality with tolerance fraction after
        # noting that all the points where the cloudmasked image
        # stopped matching after update were within 1e-16 of b7/b2 = 0.75.

        tolerance_fraction = 0.01
        checksums = [1_736_661_355, 5_997_708_807, 6_083_703_526]
        @test all(
            [abs(1 - sum(false_color_cloudmasked[i, :, :]) / checksums[i]) for i in 1:3] .< tolerance_fraction,
        )
    end

    @testset "Adaptive histogram equalization" begin
        function convert_to_255_matrix(img)::Matrix{Int}
            img_clamped = clamp.(img, 0.0, 1.0)
            return round.(Int, img_clamped * 255)
        end

        img = convert_to_255_matrix(testimage("cameraman"))
        img_eq = adapthisteq(img)
        @test sum(img_eq) == 32_387_397
    end

    @testset "Conditional adaptivehisteq" begin
        clouds = _get_false_color_cloudmasked(;
            false_color_image=false_color_image,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0,
        )

        clouds_red = clouds[1, :, :]
        clouds_red[dilated_landmask] .= 0
        tolerance_fraction = 0.01
        @test abs(1 - sum(clouds_red) / 1_320_925_065) < tolerance_fraction

        tiles = get_tiles(true_color_image; rblocks=8, cblocks=6)
        true_color_eq = conditional_histeq(true_color_image, clouds_red, tiles)

        # This differs from MATLAB script due to disparity in the implementations
        # of the adaptive histogram equalization / diffusion functions
        # For the moment testing for regression
        old_value = 6_372_159_606
        new_value = sum(to_uint8(true_color_eq[:, :, 1]))
        @test abs(1 - new_value / old_value) < 0.003

        # Use custom tile size
        side_length = size(true_color_eq, 1) ÷ 8
        tiles = get_tiles(true_color_image, side_length)
        true_color_eq = conditional_histeq(true_color_image, clouds_red, tiles)
        old_value = 6_328_796_398
        new_value = sum(to_uint8(true_color_eq[:, :, 1]))
        @test abs(1 - new_value / old_value) < 0.003
    end

    @testset "RGB to grayscale" begin
        g = rgb2gray(true_color_image)
        @test g[1000, 1000] == 92 && g[2000, 2000] == 206
    end
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
