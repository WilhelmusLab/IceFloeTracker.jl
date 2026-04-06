@testsnippet SyntheticTrackerHelpers begin
    import DataFrames: DataFrame
    import Dates: DateTime, Second
    import Images: label_components
    using IceFloeTracker

    """
        make_floe_image(n) -> BitMatrix

    Create an 8×8 BitMatrix containing a single approximately-circular floe of
    exactly `n` pixels centred in the image. Pixels are selected in order of
    increasing distance from the image centre (4.5, 4.5).
    """
    function make_floe_image(n::Int)
        center_r, center_c = 4.5, 4.5
        all_pixels = [(i, j) for i in 1:8, j in 1:8]
        sorted_pixels = sort(
            vec(all_pixels); by=p -> (p[1] - center_r)^2 + (p[2] - center_c)^2
        )
        img = falses(8, 8)
        for (r, c) in sorted_pixels[1:n]
            img[r, c] = true
        end
        return img
    end

    """
        tracker_runs_without_error(img1, time1, img2, time2) -> Bool

    Convert each BitMatrix to a labeled integer image, pass the pair to the
    FloeTracker, and return `true` if the tracker completes without error and
    produces a well-formed DataFrame.
    """
    function tracker_runs_without_error(
        img1::BitMatrix, time1::DateTime, img2::BitMatrix, time2::DateTime
    )
        labeled1 = label_components(img1)
        labeled2 = label_components(img2)
        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
            minimum_area=1,
            maximum_time_step=Second(2^19),  # 2^19 s covers the maximum test interval (2^18 s)
        )
        result = tracker([labeled1, labeled2], [time1, time2])
        return result isa DataFrame
    end

    function tracker_runs_without_error(
        floe_pixel_count::Int, dt::Int, start_time::DateTime=DateTime("2025-01-01T00:00:00")
    )
        img = make_floe_image(floe_pixel_count)
        return tracker_runs_without_error(img, start_time, img, start_time + Second(dt))
    end
end

@testitem "FloeTracker – synthetic shapes and times" setup = [SyntheticTrackerHelpers] begin
    using Dates: DateTime, Second

    base_time = DateTime("2025-01-01T00:00:00")
    time_deltas_seconds = vcat([0], [2^k for k in 0:18])

    for (floe_pixel_count, dt) in Iterators.product(1:3, time_deltas_seconds)
        @test tracker_runs_without_error(floe_pixel_count, dt, base_time) broken = true
    end

    for (floe_pixel_count, dt) in Iterators.product(4:20, time_deltas_seconds)
        @test tracker_runs_without_error(floe_pixel_count, dt, base_time)
    end
end
