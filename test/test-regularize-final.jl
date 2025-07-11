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

morph_residue = readdlm(joinpath(test_files_dir, "morph_residue.csv"), ',', Int)
local_maxima_mask = readdlm(joinpath(test_files_dir, "local_maxima_mask.csv"), ',', Int)
segment_mask = readdlm(joinpath(test_files_dir, "segment_mask.csv"), ',', Bool)
L0mask = readdlm(joinpath(test_files_dir, "L0mask.csv"), ',', Bool)
expected_regularized_holes_filled = readdlm(
    joinpath(test_files_dir, "reg_holes_filled_expected.csv"), ',', Int
)
expected_regularized_sharpened = readdlm(
    joinpath(test_files_dir, "reg_sharpened.csv"), ',', Int
)

get_final_input = readdlm(joinpath(test_files_dir, "get_final.csv"), ',', Bool)
se_erosion = se
se_dilation = se_disk2()
get_final_expected = readdlm(joinpath(test_files_dir, "get_final_expected.csv"), ',', Bool)

@testset "regularize/get_final" begin
    @testset "regularize_fill_holes/sharpening" begin
        reg_holes_filled = regularize_fill_holes(
            morph_residue, local_maxima_mask, segment_mask, L0mask, 0.3
        )

        reg_sharpened = regularize_sharpening(
            reg_holes_filled, L0mask, local_maxima_mask, segment_mask, se, 10, 2, 0.5
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

        @test expected_regularized_holes_filled == reg_holes_filled
        @test expected_regularized_sharpened == reg_sharpened
        @test expected_regularized_sharpened == reg
    end

    @testset "get_final" begin
        get_final_output = get_final(
            get_final_input, segment_mask, se_erosion, se_dilation, true
        )

        @test get_final_output == get_final_expected
    end
end