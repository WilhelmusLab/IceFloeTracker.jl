using ZipFile
using DelimitedFiles: readdlm
using IceFloeTracker: imgradientmag, to_uint8, imbinarize, adjustgamma, get_holes, se_disk4

@testset "misc. image processing" begin
    r = ZipFile.Reader("test_inputs/coins.zip")
    coins = readdlm(r.files[1], ',', Int)
    close(r)

    @testset "imgradientmag" begin
        gmag = imgradientmag(coins) |> to_uint8
        @test sum(gmag) == 2938959
    end

    @testset "imbinarize" begin
        @test sum(imbinarize(coins)) == 51638
    end

    @testset "adjustgamma" begin
        @test adjustgamma(coins, 1.5) |> sum == 5330608
    end

    @testset "get_holes" begin
        bw = coins .> 100
        @test get_holes(bw, 20, se_disk4()) |> sum == 2536
    end
end
