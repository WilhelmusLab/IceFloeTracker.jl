using IceFloeTracker: imcomplement, reconstruct, reconstruct_erosion

r = ZipFile.Reader("test_inputs/coins.zip")
coins = readdlm(r.files[1], ',', Int)
coins_gray = Gray.(coins ./ 255)
close(r)

se_disk1 = IceFloeTracker.MorphSE.StructuringElements.strel_diamond((3, 3))

_round(x) = Int(round(x, RoundNearestTiesAway))
_reconstruct(img, se, type) = reconstruct(img, se, type, false)

function run_tests(test_cases, func, se)
    for (img, expected) in test_cases
        result = func(img, se)
        @test _round(Float64(sum(result))) == _round(expected)
    end
end

@testset "reconstruct" begin
    @testset "imcomplement" begin
        @test imcomplement(coins) == 255 .- coins

        coins_gray = Gray.(coins ./ 255)
        @test imcomplement(coins_gray) == 1 .- coins_gray
    end

    @testset "open_by_reconstruction" begin
        test_cases = [(coins, 7552396), (coins_gray, 29617.239215686277)]
        run_tests(test_cases, (img, se) -> _reconstruct(img, se, "erosion"), se_disk1)
    end

    @testset "close_by_reconstruction" begin
        test_cases = [(coins, 7599858), (coins_gray, 29803.36470588235)]
        run_tests(test_cases, (img, se) -> _reconstruct(img, se, "dilation"), se_disk1)
    end

    @testset "reconstruct_erosion" begin
        test_cases = [(coins, 11179481), (coins_gray, 43841.10196078432)]
        run_tests(test_cases, reconstruct_erosion, se_disk1)
    end
end