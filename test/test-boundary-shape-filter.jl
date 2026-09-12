# Phase 6: boundary-based candidate filter.
#
# Fixtures are CLOSED (first row == last row), matching add_boundary! output.

@testsnippet BoundaryFilterSetup begin
    using DataFrames
    using IceFloeTracker.Tracking:
        BoundaryShapeDifferenceThresholdFilter,
        ChainedFilterFunction,
        PiecewiseLinearThresholdFunction,
        boundary_shape_difference,
        boundary_modified_hausdorff,
        boundary_hausdorff,
        rotate_boundary

    # closed, asymmetric L-shape
    L = [
        0.0 0.0
        3.0 0.0
        3.0 1.0
        1.0 1.0
        1.0 3.0
        0.0 3.0
        0.0 0.0
    ]
    # a deliberately dissimilar shape of comparable extent
    BAR = [
        0.0 0.0
        3.0 0.0
        3.0 0.4
        0.0 0.4
        0.0 0.0
    ]

    floe_row() = DataFrame(; boundary=[L], orientation=[0.0], area=[100.0])[1, :]
    candidate_df() = DataFrame(;
        boundary=[L, rotate_boundary(L, 0.0), BAR],
        orientation=[0.0, 0.0, 0.0],
        area=[100.0, 100.0, 100.0],
    )
end

@testitem "BoundaryShapeDifferenceThresholdFilter defaults" setup = [BoundaryFilterSetup] begin
    f = BoundaryShapeDifferenceThresholdFilter()
    @test f.metric === boundary_modified_hausdorff
    @test f.boundary_column === :boundary
    @test f.area_variable === :area
    @test f.threshold_column === :boundary_shape_difference_test
    @test f isa Function       # via AbstractFloeFilterFunction <: Function
end

@testitem "BoundaryShapeDifferenceThresholdFilter writes its score columns" setup = [BoundaryFilterSetup] begin
    f = BoundaryShapeDifferenceThresholdFilter()
    candidates = candidate_df()

    f(floe_row(), candidates, Val(:raw))

    @test "boundary_shape_difference" ∈ names(candidates)
    @test "scaled_boundary_shape_difference" ∈ names(candidates)
    @test String(f.threshold_column) ∈ names(candidates)

    # identical shapes score ~0; the dissimilar bar scores higher
    @test isapprox(candidates.boundary_shape_difference[1], 0.0; atol=1e-6)
    @test candidates.boundary_shape_difference[3] > candidates.boundary_shape_difference[1]
    @test all(isfinite, candidates.boundary_shape_difference)
end

@testitem "BoundaryShapeDifferenceThresholdFilter does not double-normalize" setup = [BoundaryFilterSetup] begin
    # The score is a distance in pixels. Unlike the mask filter there is no division by
    # :area -- count_symdiff is a raw pixel count and needs one, a Hausdorff distance
    # does not. If a division crept in, these two columns would differ by a factor of 100.
    f = BoundaryShapeDifferenceThresholdFilter()
    candidates = candidate_df()
    f(floe_row(), candidates, Val(:raw))

    @test candidates.scaled_boundary_shape_difference ==
        candidates.boundary_shape_difference
end

@testitem "BoundaryShapeDifferenceThresholdFilter subsets via the 2-arg functor" setup = [BoundaryFilterSetup] begin
    # A tight explicit threshold, so the test does not depend on the inherited default
    # bounds, which are uncalibrated for a pixel-valued score.
    f = BoundaryShapeDifferenceThresholdFilter(;
        threshold_function=PiecewiseLinearThresholdFunction(100, 700, 1e-4, 1e-4)
    )
    candidates = candidate_df()
    f(floe_row(), candidates)

    @test nrow(candidates) == 2                                  # the BAR is rejected
    @test String(f.threshold_column) ∉ names(candidates)         # test column dropped
    @test "boundary_shape_difference" ∈ names(candidates)        # score column kept
end

@testitem "BoundaryShapeDifferenceThresholdFilter honours an injected metric" setup = [BoundaryFilterSetup] begin
    a = candidate_df()
    b = candidate_df()
    BoundaryShapeDifferenceThresholdFilter(; metric=boundary_modified_hausdorff)(
        floe_row(), a, Val(:raw)
    )
    BoundaryShapeDifferenceThresholdFilter(; metric=boundary_hausdorff)(
        floe_row(), b, Val(:raw)
    )
    # max and mean nearest-neighbour distance differ for a genuinely different shape
    @test a.boundary_shape_difference[3] != b.boundary_shape_difference[3]
end

@testitem "BoundaryShapeDifferenceThresholdFilter honours boundary_column" setup = [BoundaryFilterSetup] begin
    f = BoundaryShapeDifferenceThresholdFilter(; boundary_column=:bd)
    floe = DataFrame(; bd=[L], orientation=[0.0], area=[100.0])[1, :]
    candidates = DataFrame(; bd=[L, BAR], orientation=[0.0, 0.0], area=[100.0, 100.0])

    f(floe, candidates, Val(:raw))
    @test isapprox(candidates.boundary_shape_difference[1], 0.0; atol=1e-6)
end

@testitem "BoundaryShapeDifferenceThresholdFilter composes in a chain" setup = [BoundaryFilterSetup] begin
    f = BoundaryShapeDifferenceThresholdFilter(;
        threshold_function=PiecewiseLinearThresholdFunction(100, 700, 1e-4, 1e-4)
    )
    chain = ChainedFilterFunction(; filters=[f])
    candidates = candidate_df()

    chain(floe_row(), candidates)
    @test nrow(candidates) == 2
end

# ---------------------------------------------------------------------------
# Pairing the filter with MinimumWeightMatchingFunction.
#
# :scaled_boundary_shape_difference is NOT added to the matching function's
# default `columns`. That list is validated with a subset check that, on a miss,
# logs at @debug and returns an EMPTY DataFrame (matching_functions.jl:37-42).
# So adding it there while the filter stays opt-in would make the default
# tracker silently match nothing. Opt into both together instead.
# ---------------------------------------------------------------------------

@testitem "boundary score column is absent from the matching defaults" setup = [BoundaryFilterSetup] begin
    using IceFloeTracker.Tracking: MinimumWeightMatchingFunction

    @test :scaled_boundary_shape_difference ∉ MinimumWeightMatchingFunction().columns
    @test length(MinimumWeightMatchingFunction().columns) ==
        length(MinimumWeightMatchingFunction().weights)
end

@testitem "boundary filter pairs with an opt-in matching function" setup = [BoundaryFilterSetup] begin
    using IceFloeTracker.Tracking: MinimumWeightMatchingFunction

    cols = [
        :scaled_distance,
        :relative_error_area,
        :relative_error_convex_area,
        :relative_error_major_axis_length,
        :relative_error_minor_axis_length,
        :psi_s_correlation_score,
        :scaled_shape_difference,
        :scaled_boundary_shape_difference,
    ]
    m = MinimumWeightMatchingFunction(; columns=cols, weights=ones(length(cols)))
    @test length(m.columns) == 8
    @test :scaled_boundary_shape_difference ∈ m.columns

    # the guard returns an empty frame when the score column is missing, which is
    # exactly what would happen to every default user if this column were added
    # to the defaults without also enabling the filter
    pairs = DataFrame(; head_uuid=["a"], uuid=["b"])
    @test nrow(m(pairs)) == 0
end
