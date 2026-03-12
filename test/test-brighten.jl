@testitem "brighten tests" begin
    using IceFloeTracker: get_brighten_mask, imbrighten
    import Images: Gray, AbstractGray

    @testset "get_brighten_mask" begin
        img = rand(0:255, 5, 5)
        bumped_img = img .+ 1
        mask = get_brighten_mask(img, bumped_img)
        @test all(mask .== 0)
    end

    @testset "imbrighten tests" begin
        img = [1 2; 3 4]
        brighten_mask = [1 0; 1 0]

        test_cases = [(1.25, [1 2; 4 4]), (0.1, [0 2; 0 4]), (0.9, img)]

        for (bright_factor, expected_result) in test_cases
            result = imbrighten(img, brighten_mask, bright_factor)
            @test result == expected_result
        end

        img = Gray.(img ./ maximum(img))
        result = imbrighten(img, brighten_mask .> 0, 1.25)
        # Expected result is different in this case because we aren't
        # rounding to integer precision in the end.
        expected = Gray.([1 * 1.25 2; 3 * 1.25 4] ./ 4)
        @test eltype(result) == eltype(img)
        @test result == expected
    end
end
