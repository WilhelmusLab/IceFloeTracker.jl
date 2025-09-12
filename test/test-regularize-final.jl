
@testitem "regularize/get_final" begin
    using IceFloeTracker:
        get_tiles,
        regularize_fill_holes,
        regularize_sharpening,
        _regularize,
        se_disk2,
        get_final
    using DelimitedFiles: readdlm

    se = collect(IceFloeTracker.strel_diamond((3, 3)))

    test_files_dir = joinpath(@__DIR__, "test_inputs/regularize")

    morph_residue = Gray.(readdlm(joinpath(test_files_dir, "morph_residue.csv"), ',', Int) ./ 255)
    
    local_maxima_mask = readdlm(joinpath(test_files_dir, "local_maxima_mask.csv"), ',', Int) .> 0
    segment_mask = readdlm(joinpath(test_files_dir, "segment_mask.csv"), ',', Bool)
    L0mask = readdlm(joinpath(test_files_dir, "L0mask.csv"), ',', Bool)
    
    expected_regularized_holes_filled = Gray.(readdlm(
        joinpath(test_files_dir, "reg_holes_filled_expected.csv"), ',', Int) ./ 255)

    expected_regularized_sharpened = Gray.(readdlm(
        joinpath(test_files_dir, "reg_sharpened.csv"), ',', Int) ./ 255)

    get_final_input = readdlm(joinpath(test_files_dir, "get_final.csv"), ',', Bool)
    se_erosion = se
    se_dilation = se_disk2()
    get_final_expected = readdlm(
        joinpath(test_files_dir, "get_final_expected.csv"), ',', Bool
    )

    # Differences exist only in a couple places. 
    # Potentially related to differences in the histogram equalization algorithm.

    @testset "regularize_fill_holes/sharpening" begin
        reg_holes_filled = regularize_fill_holes(
            morph_residue, local_maxima_mask, segment_mask, L0mask, 0.3
        )


        reg_sharpened = regularize_sharpening(
            expected_regularized_holes_filled, L0mask, local_maxima_mask, segment_mask, se, 10, 2, 0.5
        )

        reg = _regularize(
            morph_residue,
            local_maxima_mask,
            segment_mask,
            L0mask,
            se;
            factor=(0.3, 0.5),
            radius=10,
            amount=2,
        )

        @test maximum(abs.(expected_regularized_holes_filled .- reg_holes_filled)) .< 1/255
        # First compare as if the reg_holes_filled was exactly as expected
        @test maximum(abs.(expected_regularized_sharpened .- reg_sharpened)) .< 1/255

        # Then string both operations together and compare
        @test maximum(abs.(expected_regularized_sharpened .- reg)) .< 3/255 
    end

    @testset "get_final" begin
        get_final_output = get_final(
            get_final_input, segment_mask, se_erosion, se_dilation, true
        )

        @test get_final_output == get_final_expected
    end
end