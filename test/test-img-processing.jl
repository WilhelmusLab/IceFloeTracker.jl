using ZipFile
using DelimitedFiles: readdlm
using IceFloeTracker: imgradientmag, to_uint8, imbinarize, adjustgamma, get_holes, se_disk4

@testset "misc. image processing" begin
    r = ZipFile.Reader("test_inputs/coins.zip")
    coins = readdlm(r.files[1], ',', Int)
    close(r)

    @testset "imgradientmag" begin
        gmag = to_uint8(imgradientmag(coins))
        @test sum(gmag) == 2938959
    end

    @testset "imbinarize" begin
        @test sum(imbinarize(coins)) == 51638
    end

    @testset "adjustgamma" begin
        @test sum(adjustgamma(coins, 1.5)) == 5330608
    end

    @testset "get_holes" begin
        bw = coins .> 100
        @test sum(get_holes(bw)) == 2536
    end
end