using IceFloeTracker:
    imgradientmag,
    to_uint8,
    imbinarize,
    adjustgamma,
    get_holes,
    impose_minima,
    se_disk4,
    se_disk20,
    se_disk50

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

    @testset "impose_minima" begin
        img = readdlm("test_inputs/imposemin.csv", ',', Int)
        marker = falses(size(img))
        marker[65:70, 65:70] .= true
        @test sum(impose_minima(img, marker)) == 7675653
    end

    @testset "se_disk4" begin
        @test sum(se_disk4()) == 37
    end

    @testset "se_disk20" begin
        @test sum(se_disk20()) == 1301
    end

    @testset "se_disk50/se for landmask dilation" begin
        @test sum(se_disk50()) == 8177
    end
end
