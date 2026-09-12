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

@testitem "boundary_normalized_distance scale invariance" begin
    using IceFloeTracker.Tracking: boundary_normalized_distance

    # Two boundaries that are the same shape but different scale
    b1 = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    b2_scaled = [
        0.0 0.0
        2.0 0.0
        2.0 2.0
        0.0 2.0
    ]

    # After normalization by perimeter^2, scaled versions should have same distance
    # (or at least predictable relationship)
    dist = boundary_normalized_distance(b1, b2_scaled)

    # For identical shapes at different scales, normalized distance should be 0
    @test isapprox(dist, 0.0; atol=1e-10)
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
