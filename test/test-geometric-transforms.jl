@testitem "get_rotation_matrix identity" begin
    using IceFloeTracker.Tracking: _get_rotation_matrix

    # Zero rotation should give identity matrix
    R = _get_rotation_matrix(0.0)
    @test isapprox(R, [1.0 0.0; 0.0 1.0]; atol=1e-10)
end

@testitem "get_rotation_matrix 90 degrees" begin
    using IceFloeTracker.Tracking: _get_rotation_matrix

    # 90 degree rotation (π/2 radians)
    R = _get_rotation_matrix(π/2)
    expected = [0.0 -1.0; 1.0 0.0]
    @test isapprox(R, expected; atol=1e-10)
end

@testitem "get_rotation_matrix full circle" begin
    using IceFloeTracker.Tracking: _get_rotation_matrix

    # Full rotation should return to identity
    R = _get_rotation_matrix(2π)
    @test isapprox(R, [1.0 0.0; 0.0 1.0]; atol=1e-10)
end

@testitem "get_rotation_matrix determinant" begin
    using IceFloeTracker.Tracking: _get_rotation_matrix
    using LinearAlgebra: det

    # Rotation matrices are orthogonal with determinant 1
    for angle in [0, π/6, π/4, π/3, π/2, π, 3π/2]
        R = _get_rotation_matrix(angle)
        @test isapprox(det(R), 1.0; atol=1e-10)
    end
end

@testitem "center_boundary translation" begin
    using IceFloeTracker.Tracking: center_boundary
    using StatsBase: mean

    # Create a simple boundary (no closing point)
    boundary = [
        0.0 0.0
        2.0 0.0
        2.0 2.0
        0.0 2.0
    ]

    # Center it at origin
    centered = center_boundary(boundary; target_center=(0.0, 0.0))

    # Compute centroid of centered boundary
    centroid = vec(mean(centered; dims=1))

    # Should be approximately at origin
    @test isapprox(centroid[1], 0.0; atol=1e-6)
    @test isapprox(centroid[2], 0.0; atol=1e-6)

    # Original centroid should be at (1, 1)
    orig_centroid = vec(mean(boundary; dims=1))
    @test isapprox(orig_centroid[1], 1.0; atol=1e-6)
    @test isapprox(orig_centroid[2], 1.0; atol=1e-6)
end

@testitem "center_boundary to arbitrary target" begin
    using IceFloeTracker.Tracking: center_boundary
    using StatsBase: mean

    boundary = [
        0.0 0.0
        1.0 0.0
        1.0 1.0
        0.0 1.0
    ]

    target = (5.0, 10.0)
    centered = center_boundary(boundary; target_center=target)

    centroid = vec(mean(centered; dims=1))
    @test isapprox(centroid[1], target[1]; atol=1e-6)
    @test isapprox(centroid[2], target[2]; atol=1e-6)
end

@testitem "rotate_boundary preserves boundary shape" begin
    using IceFloeTracker.Tracking: rotate_boundary
    using StatsBase: mean

    # Square boundary (no closing point)
    boundary = [
        1.0 0.0
        0.0 1.0
        -1.0 0.0
        0.0 -1.0
    ]

    # Rotate by 90 degrees (π/2)
    rotated = rotate_boundary(boundary, π/2)

    # Shape should be preserved (same number of points)
    @test size(rotated, 1) == size(boundary, 1)

    # Centroid should stay the same (rotation around centroid)
    orig_centroid = vec(mean(boundary; dims=1))
    rot_centroid = vec(mean(rotated; dims=1))
    @test isapprox(orig_centroid[1], rot_centroid[1]; atol=1e-6)
    @test isapprox(orig_centroid[2], rot_centroid[2]; atol=1e-6)
end

@testitem "rotate_boundary with explicit center" begin
    using IceFloeTracker.Tracking: rotate_boundary

    # Simple 2-point boundary
    boundary = [
        1.0 0.0
        2.0 0.0
    ]

    # Rotate around (0, 0) by 90 degrees
    center = (0.0, 0.0)
    rotated = rotate_boundary(boundary, π/2; center=center)

    # First point (1, 0) should rotate to (0, 1)
    @test isapprox(rotated[1, 1], 0.0; atol=1e-6)
    @test isapprox(rotated[1, 2], 1.0; atol=1e-6)
end

@testitem "rotate_boundary 360 degrees returns original" begin
    using IceFloeTracker.Tracking: rotate_boundary

    boundary = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]

    # Full rotation
    full_rotation = rotate_boundary(boundary, 2π)

    # Should be approximately the same
    for i in axes(boundary, 1)
        @test isapprox(full_rotation[i, 1], boundary[i, 1]; atol=1e-5)
        @test isapprox(full_rotation[i, 2], boundary[i, 2]; atol=1e-5)
    end
end
