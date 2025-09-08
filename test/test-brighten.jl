
@testitem "brighten tests" begin
    using IceFloeTracker: get_brighten_mask, imbrighten

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
    end
end