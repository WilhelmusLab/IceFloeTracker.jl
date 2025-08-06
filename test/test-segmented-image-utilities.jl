@ntestset "$(@__FILE__)" begin
    @ntestset "segmentation_summary" begin
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
    @ntestset "segmentation_comparison" begin
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
        
        @ntestset "Self-similar results with zero labels" begin
            for label_map in [
                [0],
                [0 0],
                [0 0 0],
                [
                    0 0 0
                    0 0 0
                    0 0 0
                ],
            ]
                @show label_map
                test_identical_segmentation_properties(
                    label_map;
                    normalized_validated_area=0.0,
                    normalized_measured_area=0.0,
                    fractional_intersection=NaN,
                )
            end
        end
        @ntestset "Self-similar results with all non-zero labels" begin
            for label_map in [[1], [1, 2], [1 2 3]]
                @show label_map
                test_identical_segmentation_properties(
                    label_map;
                    normalized_validated_area=1,
                    normalized_measured_area=1,
                    fractional_intersection=1,
                )
            end
        end
        @ntestset "Different label indices" begin
            # The checks should be agnostic to what the label indices are
            for (label_map_1, label_map_2) in
                [([1], [2]), ([1, 2], [2, 1]), ([1 2 3], [3 1 2])]
                @show label_map_1, label_map_2
                test_segmentation_properties(
                    label_map_1,
                    label_map_2;
                    normalized_validated_area=1,
                    normalized_measured_area=1,
                    fractional_intersection=1,
                )
            end
        end
        @ntestset "Overlay fractions" begin
            test_segmentation_properties([0 0], [0 1]; fractional_intersection=NaN)
            test_segmentation_properties([0 1], [0 0]; fractional_intersection=0.0)
            test_segmentation_properties([0 1], [1 0]; fractional_intersection=0.0)
            test_segmentation_properties([0 1], [0 1]; fractional_intersection=1.0)
            test_segmentation_properties([0 1], [1 1]; fractional_intersection=1.0)
            test_segmentation_properties([1 1], [0 1]; fractional_intersection=0.5)
            test_segmentation_properties([1 1], [1 1]; fractional_intersection=1.0)
            test_segmentation_properties([1 1 1], [1 0 0]; fractional_intersection=1 / 3)
            test_segmentation_properties([1 1 1], [1 1 0]; fractional_intersection=2 / 3)
            test_segmentation_properties([1 1 1], [1 1 1]; fractional_intersection=3 / 3)
            test_segmentation_properties([1 0 0], [1 1 1]; fractional_intersection=1)
            test_segmentation_properties([1 1 0], [1 1 1]; fractional_intersection=1)
            test_segmentation_properties([1 1 1], [1 1 1]; fractional_intersection=1)
        end
    end
end
