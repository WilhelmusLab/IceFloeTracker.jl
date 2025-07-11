using IceFloeTracker:
    _get_false_color_cloudmasked,
    convert_to_255_matrix,
    adapthisteq,
    conditional_histeq,
    rgb2gray,
    to_uint8,
    histeq

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
    dilated_landmask = BitMatrix(load(joinpath(datadir, "matlab_landmask.png")))
end

function test_cloud_image_workflow()
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
        @test all([abs(1 - sum(false_color_cloudmasked[i, :, :])/checksums[i]) for i in 1:3] .< tolerance_fraction)
    end
end

function test_adaphisteq()
    @testset "Adaptive histogram equalization" begin
        img = convert_to_255_matrix(testimage("cameraman"))
        img_eq = adapthisteq(img)
        @test sum(img_eq) == 32_387_397
    end
end

function test_conditional_adaptivehisteq()
    @testset "Conditional adaptivehisteq" begin
        clouds = _get_false_color_cloudmasked(;
            false_color_image=false_color_image,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0
        )

        clouds_red = clouds[1, :, :]
        clouds_red[dilated_landmask] .= 0
        tolerance_fraction = 0.01
        @test abs(1 - sum(clouds_red) / 1_320_925_065) < tolerance_fraction

        # Using rblocks = 8, cblocks = 6
        true_color_eq = conditional_histeq(true_color_image, clouds_red, 8, 6)

        # This differs from MATLAB script due to disparity in the implementations
        # of the adaptive histogram equalization / diffusion functions
        # For the moment testing for regression
        @test sum(to_uint8(true_color_eq[:, :, 1])) == 6_372_159_606

        # Use custom tile size
        side_length = size(true_color_eq, 1) ÷ 8
        true_color_eq = conditional_histeq(true_color_image, clouds_red, side_length)
        @test sum(to_uint8(true_color_eq[:, :, 1])) == 6_328_796_398
    end
end

function test_rgb2gray()
    @testset "RGB to grayscale" begin
        g = rgb2gray(true_color_image)
        @test g[1000, 1000] == 92 && g[2000, 2000] == 206
    end
end

function test_histeq()
    @testset "histeq" begin
        imgs = [
            [
                4 4 4 4 4
                3 4 5 4 3
                3 5 5 5 3
                3 4 5 4 3
                4 4 4 4 4
            ],
            # Edge cases
            # maximum a power of 2
            [
                0 0 0
                2 2 2
            ],
            # maximum at 255
            [
                0 0 0
                255 255 255
            ],
        ]

        _exp = [
            128 128 128
            255 255 255
        ]

        expected = [
            [
                204 204 204 204 204
                61 204 255 204 61
                61 255 255 255 61
                61 204 255 204 61
                204 204 204 204 204
            ],
            _exp,
            _exp,
        ]

        @test all(histeq(imgs[i]) == expected[i] for i in 1:3)
    end
end

@testset "Conditional adaptivehisteq" begin
    test_cloud_image_workflow()
    test_adaphisteq()
    test_conditional_adaptivehisteq()
    test_rgb2gray()
    test_histeq()
end
