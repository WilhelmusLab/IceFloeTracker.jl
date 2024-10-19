# WIP
using ZipFile # delete later
using DelimitedFiles: readdlm
using IceFloeTracker: imcomplement

r = ZipFile.Reader("test/test_inputs/coins.zip")
coins = readdlm(r.files[1], ',', Int)
close(r)

@testset "reconstruct" begin
    @testset "imcomplement" begin
        @test imcomplement(coins) == 255 .- coins

        coins_gray = Gray.(coins ./ 255)
        @test imcomplement(coins_gray) == 1 .- coins_gray
    end

    @testset "open_by_reconstruction" begin
        se_disk1 = IceFloeTracker.MorphSE.StructuringElements.strel_diamond((3, 3))
        coins_opened = open_by_reconstruction(coins, se_disk1)
        @test sum(coins_opened) == 7552396
    end

    @testset "reconstruct_erosion" begin
        coins_r_erosion = reconstruct_erosion(coins, se_disk1)
        @test sum(coins_r_erosion) == 11179481
    end
end
