@testitem "MinimumWeightMatchingFunction happy path" begin
    using IceFloeTracker.Tracking: MinimumWeightMatchingFunction
    using DataFrames

    minimum_weight_matching_function = MinimumWeightMatchingFunction(; columns=[:dx])

    # Happy path for a single head_uuid with a clear minimum weight match.
    @test minimum_weight_matching_function(
        DataFrame(; head_uuid=["a", "a"], uuid=["c", "d"], dx=[1.0, 2.0])
    ) == DataFrame(; head_uuid=["a"], uuid=["c"], dx=[1.0], w=[1.0])

    @test minimum_weight_matching_function(
        DataFrame(;
            head_uuid=["a", "a", "b", "b"],
            uuid=["c", "d", "e", "f"],
            dx=[1.0, 2.0, 3.0, 4.0],
        ),
    ) == DataFrame(; head_uuid=["a", "b"], uuid=["c", "e"], dx=[1.0, 3.0], w=[1.0, 3.0])

    # If there is a tie, only one match is returned for each head_uuid, and it should be the one with the lowest dx.
    @test minimum_weight_matching_function(
        DataFrame(;
            head_uuid=["a", "a", "b", "b"],
            uuid=["c", "d", "c", "d"],
            dx=[1.0, 2.0, 3.0, 4.0],
        ),
    ) == DataFrame(; head_uuid=["a"], uuid=["c"], dx=[1.0], w=[1.0])
    @test minimum_weight_matching_function(
        DataFrame(;
            head_uuid=["a", "a", "b", "b"],
            uuid=["c", "d", "c", "d"],
            dx=[4.0, 3.0, 2.0, 1.0],
        ),
    ) == DataFrame(; head_uuid=["b"], uuid=["d"], dx=[1.0], w=[1.0])

    # Ideally, the matching will be clearly defined
    @test minimum_weight_matching_function(
        DataFrame(;
            head_uuid=["a", "a", "b", "b"],
            uuid=["c", "d", "c", "d"],
            dx=[1.0, 2.0, 2.0, 1.0],
        ),
    ) == DataFrame(; head_uuid=["a", "b"], uuid=["c", "d"], dx=[1.0, 1.0], w=[1.0, 1.0])
end
@testitem "MinimumWeightMatchingFunction empty DataFrame" begin
    using IceFloeTracker.Tracking: MinimumWeightMatchingFunction
    using DataFrames

    minimum_weight_matching_function = MinimumWeightMatchingFunction(; columns=[:dx])

    # Test that MinimumWeightMatchingFunction works on an empty DataFrame with the correct columns. 
    empty_dataframe = DataFrame(
        map(
            col -> col => String[],
            [minimum_weight_matching_function.columns..., :head_uuid, :uuid],
        ),
    )
    minimum_weight_matching_function(empty_dataframe)
end