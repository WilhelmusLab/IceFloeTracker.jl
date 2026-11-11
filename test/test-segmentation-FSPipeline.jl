@testitem "FSPipeline.Segment – simple case" tags = [:e2e, :smoke] begin
    import DataFrames: DataFrame, nrow
    import Images: RGB
    dataset = Watkins2026Dataset(; ref="v0.2")

    case = first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset))
    segments = FSPipeline.Segment()(
        RGB.(modis_truecolor(case)), RGB.(modis_falsecolor(case)), modis_landmask(case)
    )
    expected_segment_count = validated_floe_properties(case) |> DataFrame |> nrow
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 0.7
    # check that the larger-than-image tile size case works
    segments = FSPipeline.Segment(tile_size_pixels=500)(
        RGB.(modis_truecolor(case)), RGB.(modis_falsecolor(case)), modis_landmask(case)
    )
    @test length(segments.segment_labels) ≈ expected_segment_count rtol = 0.7
end

@testitem "FSPipeline.Segment – sample of cases" setup = [Segmentation] tags = [:e2e] begin
    import StatsBase: mean
    dataset = Watkins2026Dataset(; ref="v0.2")
    results = run_and_validate_segmentation(
        filter(
            c -> (
                c.visible_floes == "yes" &&
                c.cloud_fraction_manual < 0.5 &&
                c.case_number % 3 == 0
            ),
            dataset,
        ),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )
    @test all(results.success)
    save("./test_outputs/sample_error_metrics_FSPipeline.csv", results)
    # Aggregate performance measures
    mean_recall = round(mean(skipnanormissing(results.recall)); digits=2)
    mean_precision = round(mean(skipnanormissing(results.precision)); digits=2)
    mean_F_score = round(mean(skipnanormissing(results.F_score)); digits=2)

    # Good performance might look like this:
    @test mean_recall ≥ 0.9 broken = true
    @test mean_precision ≥ 0.9 broken = true
    @test mean_F_score ≥ 0.9 broken = true

    # Better performance might look like this:
    @test mean_recall ≥ 0.8 broken = true # Note: Increase in recall without increase in precision may indicate higher false-positive rate.
    @test mean_precision ≥ 0.8 broken = true
    @test mean_F_score ≥ 0.8 broken = true

    # Current performance should look at least as good as this:
    @test mean_recall ≥ 0.6
    @test mean_precision ≥ 0.3
    @test round(mean_F_score; digits=1) ≥ 0.4

    # return current performance
    @show mean_recall
    @show mean_precision
    @show mean_F_score
end

@testitem "FSPipeline.Segment – detailed tests" setup = [Segmentation] tags = [:e2e] begin
    dataset = Watkins2026Dataset(; ref="v0.2")
    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 6 && c.satellite == "terra"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )

    @test 0.36 ≈ labeled_fraction atol = 0.1
    @test 0.88 ≤ round(recall; digits=2)
    @test 0.72 ≤ round(precision; digits=2)
    @test 0.79 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 14 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )

    @test 0.16 ≈ labeled_fraction atol = 0.1
    @test 0.79 ≤ round(recall; digits=2)
    @test 0.59 ≤ round(precision; digits=2)
    @test 0.7 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 61 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )

    @test 0.23 ≈ labeled_fraction atol = 0.1
    @test 0.36 ≤ round(recall; digits=2)
    @test 0.84 ≤ round(precision; digits=2)
    @test 0.53 ≤ round(F_score; digits=2)

    (; labeled_fraction, recall, precision, F_score) = run_and_validate_segmentation(
        first(filter(c -> (c.case_number == 63 && c.satellite == "aqua"), dataset)),
        FSPipeline.Segment();
        output_directory="./test_outputs/",
    )
    # Note: Validation dataset currently doesn't include the floes intersecting the edge.
    # Improving the segmentation lowered the scores here due to these floes.

    @test labeled_fraction ≈ 0.64 rtol = 0.1
    @test 0.90 ≤ round(recall; digits=2)
    @test 0.91 ≤ round(precision; digits=2)
    @test 0.93 ≤ round(F_score; digits=2)
end

@testitem "FSPipeline.Segment – image types" setup = [Segmentation] tags = [:e2e] begin
    import Images: RGB, RGBA, n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64
    dataset = Watkins2026Dataset(; ref="v0.2")
    case::Case = first(filter(c -> (c.case_number == 6 && c.satellite == "aqua"), dataset))
    algorithm = FSPipeline.Segment()
    baseline = run_and_validate_segmentation(
        case, algorithm; output_directory="./test_outputs/"
    )

    paired_types = [n0f8, n6f10, n4f12, n2f14, n0f16, float32, float64]
    unary_types = [RGB, RGBA, paired_types...]

    for T in unary_types
        @test results_invariant_for(T; baseline, algorithm, case)
    end

    for T in paired_types
        @test results_invariant_for(RGB, T; baseline, algorithm, case)
    end
end

@testitem "FSPipeline.objectwise_compare_segmentation: reflectance columns" begin
    import IceFloeTracker: expand_labels
    import DataFrames: nrow
    import Images: Gray
    import StatsBase: mean

    objectwise_compare_segmentation = FSPipeline.objectwise_compare_segmentation

    # Two segmentations of the same two floes; floe 1 is one row taller in the second.
    indexmap1 = zeros(Int, 12, 12)
    indexmap1[2:4, 2:4] .= 1
    indexmap1[7:10, 7:10] .= 2
    indexmap2 = copy(indexmap1)
    indexmap2[5, 2:4] .= 1
    img = Gray.(reshape(range(0, 1; length=144), 12, 12))
    radius = 2

    # Mean over each region and over its external ring, 0 when the ring is empty.
    function reference_means(indexmap, labels)
        ring = expand_labels(indexmap, radius) .- indexmap
        region_mean = [mean(img[indexmap .== L]) for L in labels]
        ring_mean = [any(ring .== L) ? mean(img[ring .== L]) : 0 for L in labels]
        return region_mean, ring_mean
    end

    df = objectwise_compare_segmentation(indexmap1, indexmap2, img; expand_radius=radius)
    @test nrow(df) == 2
    for (prefix, indexmap) in (("s1", indexmap1), ("s2", indexmap2))
        region_mean, ring_mean = reference_means(indexmap, df[:, "$(prefix)_label"])
        @test all(df[:, "$(prefix)_reflectance_mean"] .≈ region_mean)
        @test all(df[:, "$(prefix)_reflectance_bdry_mean"] .≈ ring_mean)
        @test all(df[:, "$(prefix)_reflectance_bdry_contrast"] .≈ region_mean .- ring_mean)
    end

    # Floes covering the whole image have an empty ring: the boundary mean is 0.
    full = [1 1 2 2; 1 1 2 2]
    full_img = Gray.([0.1 0.2 0.3 0.4; 0.5 0.6 0.7 0.8])
    df = objectwise_compare_segmentation(full, full, full_img; expand_radius=radius)
    @test nrow(df) == 2
    @test all(iszero, df.s1_reflectance_bdry_mean)
    @test all(iszero, df.s2_reflectance_bdry_mean)
end
