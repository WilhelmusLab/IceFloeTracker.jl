@testitem "segmentation_summary" begin
    using Images: SegmentedImage
    """
    Test the segmentation measurements for one label map.
    """
    function test_segmentation_summary_properties(label_map; kwargs...)
        segmented = SegmentedImage(label_map, label_map)
        summary = segmentation_summary(segmented)
        for (key, value) in kwargs
            @test isequal(summary[key], value)
        end
    end

    test_segmentation_summary_properties([0]; labeled_fraction=0.0)
    test_segmentation_summary_properties([1]; labeled_fraction=1.0)
    test_segmentation_summary_properties([0 0]; labeled_fraction=0.0)
    test_segmentation_summary_properties([0 1]; labeled_fraction=0.5)
    test_segmentation_summary_properties([1 0]; labeled_fraction=0.5)
    test_segmentation_summary_properties([1 1]; labeled_fraction=1.0)
    test_segmentation_summary_properties([1 0 0]; labeled_fraction=1 / 3)
    test_segmentation_summary_properties([1 1 0]; labeled_fraction=2 / 3)
    test_segmentation_summary_properties(
        [
            1 1
            1 1
        ];
        labeled_fraction=1.0,
    )
end
@testitem "segmentation_comparison" begin
    using Images: SegmentedImage

    """
    Test the segmentation comparison for two label maps.
    """
    function test_segmentation_properties(label_map_1, label_map_2; kwargs...)
        segmented_1 = SegmentedImage(label_map_1, label_map_1)
        segmented_2 = SegmentedImage(label_map_2, label_map_2)
        comparison = segmentation_comparison(segmented_1, segmented_2)
        for (key, value) in kwargs
            @test isequal(comparison[key], value)
        end
    end

    """
    Test the segmentation comparison for the same label_map.
    """
    function test_identical_segmentation_properties(label_map; kwargs...)
        return test_segmentation_properties(label_map, label_map; kwargs...)
    end

    # Self-similar results with zero labels
    test_identical_segmentation_properties([0]; recall=NaN)
    test_identical_segmentation_properties([0 0]; recall=NaN)
    test_identical_segmentation_properties([0 0 0]; recall=NaN)
    test_identical_segmentation_properties(
        [
            0 0 0
            0 0 0
            0 0 0
        ];
        recall=NaN,
    )

    # Self-similar results with all non-zero labels
    test_identical_segmentation_properties([1]; recall=1)
    test_identical_segmentation_properties([1 2]; recall=1)
    test_identical_segmentation_properties([1 2 3]; recall=1)

    # Different label indices
    # The checks should be agnostic to what the label indices are
    test_segmentation_properties([1], [2]; recall=1)
    test_segmentation_properties([1, 2], [2, 1]; recall=1)
    test_segmentation_properties([1 2 3], [3 1 2]; recall=1)

    # recall, precision and F_score
    test_segmentation_properties([0 0], [0 1]; recall=NaN, precision=0, F_score=NaN)
    test_segmentation_properties([0 1], [0 0]; recall=0.0, precision=NaN, F_score=NaN)
    test_segmentation_properties([0 1], [1 0]; recall=0.0, precision=0.0, F_score=NaN)
    test_segmentation_properties([0 1], [0 1]; recall=1.0, precision=1.0, F_score=1.0)
    test_segmentation_properties(
        [0 1], [1 1]; recall=1.0, precision=0.5, F_score=2 * (1 * 0.5) / (1 + 0.5)
    )
    test_segmentation_properties(
        [1 1], [0 1]; recall=0.5, precision=1.0, F_score=2 * (1 * 0.5) / (1 + 0.5)
    )
    test_segmentation_properties([1 1], [1 1]; recall=1.0, precision=1.0)
    test_segmentation_properties(
        [1 1 1],
        [1 0 0];
        recall=1 / 3,
        precision=1.0,
        F_score=2 * (1 * (1 / 3)) / (1 + (1 / 3)),
    )
    test_segmentation_properties([1 1 1], [1 1 0]; recall=2 / 3, precision=1.0)
    test_segmentation_properties([1 1 1], [1 1 1]; recall=3 / 3, precision=1.0)
    test_segmentation_properties([1 0 0], [1 1 1]; recall=1, precision=1 / 3)
    test_segmentation_properties([1 1 0], [1 1 1]; recall=1, precision=2 / 3)
    test_segmentation_properties([1 1 1], [1 1 1]; recall=1, precision=1.0)
end

@testitem "Stitch clusters" begin

    # test stitch_clusters by creating an image indexmap with
    # a rectangle divided into 4 clusters, and use the stitch_clusters
    # function to put it back together again

    import IceFloeTracker: get_tiles, stitch_clusters
    import Images: SegmentedImage

    test_im = zeros(Int64, (10, 10))
    test_im[2:5, 2:5] .= 1
    test_im[6:9, 2:5] .= 2
    test_im[2:5, 6:9] .= 3
    test_im[6:9, 6:9] .= 4
    tiles = get_tiles(test_im, 5) # divide into 4 tiles
    segments = SegmentedImage(ones(size(test_im)), test_im)
    stitched_segments = stitch_clusters(segments, tiles)
    @test all(stitched_segments[2:9, 2:9] .== 1)
end
