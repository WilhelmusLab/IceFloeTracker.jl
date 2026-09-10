@testitem "boundary_mse_aligned identity" begin
    using IceFloeTracker.Tracking: boundary_mse_aligned

    # Distance from boundary to itself should be zero
    boundary = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    dist = boundary_mse_aligned(boundary, boundary)
    @test isapprox(dist, 0.0; atol=1e-10)
end

@testitem "boundary_mse_aligned symmetry" begin
    using IceFloeTracker.Tracking: boundary_mse_aligned

    b1 = [
        0.0 0.0
        2.0 0.0
        2.0 2.0
        0.0 2.0
    ]

    b2 = [
        1.0 1.0
        3.0 1.0
        3.0 3.0
        1.0 3.0
    ]

    dist12 = boundary_mse_aligned(b1, b2)
    dist21 = boundary_mse_aligned(b2, b1)

    @test isapprox(dist12, dist21; atol=1e-10)
end

@testitem "boundary_mse_aligned translation invariance" begin
    using IceFloeTracker.Tracking: boundary_mse_aligned, center_boundary

    boundary = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    # Create a translated copy
    translated = boundary .+ [5.0 7.0]

    # Distance should be same regardless of translation (both are translated together)
    # Actually, after centering, they should be identical
    # So distance should be zero
    dist = boundary_mse_aligned(boundary, translated)
    @test isapprox(dist, 0.0; atol=1e-10)
end

@testitem "boundary_mse_aligned increases with difference" begin
    using IceFloeTracker.Tracking: boundary_mse_aligned

    b1 = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    b2_small_diff = [
        0.1 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    b2_large_diff = [
        0.5 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    dist_small = boundary_mse_aligned(b1, b2_small_diff)
    dist_large = boundary_mse_aligned(b1, b2_large_diff)

    @test dist_small < dist_large
end

@testitem "boundary_normalized_distance identity" begin
    using IceFloeTracker.Tracking: boundary_normalized_distance

    boundary = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    dist = boundary_normalized_distance(boundary, boundary)
    @test isapprox(dist, 0.0; atol=1e-10)
end

@testitem "boundary_normalized_distance symmetry" begin
    using IceFloeTracker.Tracking: boundary_normalized_distance

    b1 = [
        0.0 0.0
        2.0 0.0
        2.0 2.0
        0.0 2.0
    ]

    b2 = [
        1.0 1.0
        3.0 1.0
        3.0 3.0
        1.0 3.0
    ]

    dist12 = boundary_normalized_distance(b1, b2)
    dist21 = boundary_normalized_distance(b2, b1)

    @test isapprox(dist12, dist21; atol=1e-10)
end

@testitem "boundary_normalized_distance translation invariance" begin
    using IceFloeTracker.Tracking: boundary_normalized_distance

    boundary = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    translated = boundary .+ [5.0 7.0]
    dist = boundary_normalized_distance(boundary, translated)
    @test isapprox(dist, 0.0; atol=1e-10)
end

@testitem "boundary_normalized_distance is dimensionless, not scale-invariant" begin
    using IceFloeTracker.Tracking: boundary_normalized_distance, boundary_mse_aligned

    unit = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]
    doubled = 2 .* unit

    # Dividing by mean perimeter squared makes the score dimensionless so that a
    # single threshold can span floe sizes -- it does NOT make it scale-invariant.
    # The same shape at two scales is a real difference, which for floe tracking is
    # signal (a floe that doubled in size is probably not the same floe), not noise.
    dist = boundary_normalized_distance(unit, doubled)
    @test dist > 0
    @test isapprox(dist, 0.0247; atol=1e-4)

    # and it is genuinely a normalization: the raw MSE is far larger
    @test boundary_mse_aligned(unit, doubled) > dist

    # scaling BOTH inputs by the same factor leaves the score unchanged, which is
    # the invariance this metric actually provides
    @test isapprox(
        boundary_normalized_distance(unit, doubled),
        boundary_normalized_distance(10 .* unit, 10 .* doubled);
        atol=1e-12,
    )
end

@testitem "boundary_euclidean_distance identity" begin
    using IceFloeTracker.Tracking: boundary_euclidean_distance

    boundary = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    dist = boundary_euclidean_distance(boundary, boundary)
    @test isapprox(dist, 0.0; atol=1e-10)
end

@testitem "boundary_euclidean_distance symmetry" begin
    using IceFloeTracker.Tracking: boundary_euclidean_distance

    b1 = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
    ]

    b2 = [
        0.1 0.1
        1.1 0.0
        1.1 1.1
    ]

    dist12 = boundary_euclidean_distance(b1, b2)
    dist21 = boundary_euclidean_distance(b2, b1)

    @test isapprox(dist12, dist21; atol=1e-10)
end

@testitem "boundary_euclidean_distance requires same length" begin
    using IceFloeTracker.Tracking: boundary_euclidean_distance

    b1 = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
    ]

    b2 = [
        0.0 0.0
        1.0 0.0
    ]

    # Should error or handle gracefully when lengths don't match
    # (test implementation will determine behavior)
    @test_throws Exception boundary_euclidean_distance(b1, b2)
end

@testitem "boundary_euclidean_distance point-wise distance" begin
    using IceFloeTracker.Tracking: boundary_euclidean_distance

    # Simple test: distance between shifted points
    b1 = [
        0.0 0.0
        1.0 0.0
    ]

    b2 = [
        0.0 1.0
        1.0 1.0
    ]

    # Each point is 1 unit away, 2 points, should sum to 2
    dist = boundary_euclidean_distance(b1, b2)
    @test isapprox(dist, 2.0; atol=1e-10)
end
