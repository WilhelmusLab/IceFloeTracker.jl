using IceFloeTracker: get_tiles, regularize_fill_holes, regularize_sharpening, _regularize
using DelimitedFiles: readdlm

se = collect(IceFloeTracker.MorphSE.strel_diamond((3, 3)))

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
