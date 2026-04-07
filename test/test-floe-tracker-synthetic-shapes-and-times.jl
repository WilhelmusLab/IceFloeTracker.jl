@testitem "FloeTracker – synthetic shapes" setup = [TrackerValidation] begin
    using Dates: DateTime
    function tracker_runs_without_error(
        tracker::AbstractTracker,
        img::Matrix{Int};
        start=DateTime("2025-01-01T00:00:00"),
        end_=DateTime("2025-01-01T00:00:01"),
    )
        return is_wellformed_tracker_result(tracker([img, img], [start, end_]))
    end

    tracker = FloeTracker(;
        filter_function=FilterFunction(),
        matching_function=MinimumWeightMatchingFunction(),
        minimum_area=1,
    )

    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0
            0 0 0
            0 0 0
        ],
    )

    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0
            0 1 0
            0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/911

    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0
            0 1 0
            0 1 0
            0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 0 0
            0 0 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 0 1 0
            0 1 0 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 1 0
            0 0 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 0 0
            0 1 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 0 1 0 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 1 0
            0 1 1 0
            0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 0 0
            0 1 1 0
            0 1 0 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 1 1 1 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 0 1 0 0
            0 1 1 1 1 0
            0 0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0 0
            0 1 1 1 1 1 0
            0 0 0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 0 0
            0 1 0 0 0
            0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 1 1 1 0
            0 1 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 0
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 0 0
            0 0 1 0 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 0 0
            0 1 1 1 0
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 0 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 0 1 0 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 0 1 1 1 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 1
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 1 1 1
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 1 1 1 0 0
            0 0 1 0 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 0 0
            0 1 1 1 0 0
            0 1 1 1 1 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 1 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 1 1 1 0 0
            0 0 0 0 0 0
        ],
    )
end