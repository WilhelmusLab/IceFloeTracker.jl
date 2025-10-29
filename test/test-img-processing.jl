
@testitem "misc. image processing" begin
    using IceFloeTracker.LopezAcosta2019Tiling:
        imgradientmag, adjustgamma, get_holes, impose_minima
    using Images
    using ZipFile
    import DelimitedFiles: readdlm

    r = ZipFile.Reader("test_inputs/coins.zip")
    coins = readdlm(r.files[1], ',', Int)
    close(r)

    @testset "imgradientmag" begin
        gmag = to_uint8(imgradientmag(coins))
        @test sum(gmag) == 2938959
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

    @testset "impose_minima" begin
        img = readdlm("test_inputs/imposemin.csv", ',', Int)
        marker = falses(size(img))
        marker[65:70, 65:70] .= true
        @test sum(impose_minima(img, marker)) == 7675653
    end
end
