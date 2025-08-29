
@testitem "reconstruct" begin
    using IceFloeTracker: imcomplement, reconstruct, strel_diamond
    using ZipFile

    r = ZipFile.Reader("test_inputs/coins.zip")
    coins = readdlm(r.files[1], ',', Int)
    close(r)

    se_disk1 = strel_diamond((3, 3))

    _round(x) = Int(round(x, RoundNearestTiesAway))
    _reconstruct(img, se, type) = reconstruct(img, se, type, false)

    function run_tests(test_cases, func, se)
        for (img, expected) in test_cases
            result = func(img, se)
            @test _round(Float64(sum(result))) == _round(expected)
        end
    end

    @testset "imcomplement" begin
        @test imcomplement(coins) == 255 .- coins

        coins_gray = Gray.(coins ./ 255)
        @test imcomplement(coins_gray) == 1 .- coins_gray
    end

    @testset "open_by_reconstruction" begin
        test_cases = [(coins, 7552396)]
        run_tests(test_cases, (img, se) -> _reconstruct(img, se, "erosion"), se_disk1)
    end

    @testset "close_by_reconstruction" begin
        test_cases = [(coins, 7599858)]
        run_tests(test_cases, (img, se) -> _reconstruct(img, se, "dilation"), se_disk1)
    end
end
