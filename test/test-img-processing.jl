
@testitem "misc. image processing" begin
    using IceFloeTracker.LopezAcosta2019Tiling: adjustgamma, get_holes
    using IceFloeTracker.ImageUtils: binarize_mask
    using Images
    using ZipFile
    import DelimitedFiles: readdlm

    r = ZipFile.Reader("test_inputs/coins.zip")
    coins = readdlm(r.files[1], ',', Int)
    close(r)

    @testset "imgradientmag" begin
        gmag = to_uint8(imgradientmag(coins))
        @test sum(gmag) == 442528
    end

    @testset "image binarization" begin
        f = AdaptiveThreshold(coins)
        @test sum(binarize(coins, f)) == 51638
    end

    @testset "adjustgamma" begin
        @test sum(adjustgamma(coins, 1.5)) == 5330608
    end

    @testset "get_holes" begin
        bw = coins .> 100
        @test sum(get_holes(bw)) == 2536
    end

    @testset "binarize_mask" begin
        low_nonzero = Gray.(N0f8.([0.0 1 / 255; 2 / 255 0.0]))

        @test binarize_mask(low_nonzero) == BitMatrix([0 1; 1 0])
        @test binarize_mask(low_nonzero; tol=0.01) == BitMatrix([0 0; 0 0])
        @test binarize_mask(low_nonzero; tol=1 / 255) == BitMatrix([0 0; 1 0])
    end
end
