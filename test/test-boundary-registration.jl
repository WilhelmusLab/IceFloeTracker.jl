# Phase 4: boundary-curve analogues of shape_difference_rotation / register.
#
# Fixtures here are CLOSED (first point repeated as last), matching what
# add_boundary! actually produces. They are also ASYMMETRIC under rotation --
# a square would make the recovered angle ambiguous at 90 degree multiples.

@testsnippet BoundaryRegSetup begin
    using IceFloeTracker.Tracking:
        shape_difference_rotation_boundary,
        register_boundary,
        rotate_boundary,
        boundary_normalized_distance,
        boundary_mse_aligned,
        register_default_angles_rad,
        prior_test_angles

    # closed L-shape: no rotational symmetry, so the optimum angle is unique
    L = [
        0.0 0.0
        3.0 0.0
        3.0 1.0
        1.0 1.0
        1.0 3.0
        0.0 3.0
        0.0 0.0
    ]
end

@testitem "shape_difference_rotation_boundary return shape" setup = [BoundaryRegSetup] begin
    angles = [0.0, 0.1, 0.2, 0.3]
    result = shape_difference_rotation_boundary(L, L, angles)

    @test length(result) == length(angles)
    # mirrors the mask-based shape_difference_rotation contract
    @test all(r -> propertynames(r) == (:angle, :shape_difference), result)
    # angles come back in the order supplied, unmodified
    @test [r.angle for r in result] == angles
    @test all(r -> r.shape_difference isa Real, result)
end

@testitem "shape_difference_rotation_boundary identity is minimal at zero" setup = [BoundaryRegSetup] begin
    angles = [deg2rad(d) for d in -20:5:20]
    result = shape_difference_rotation_boundary(L, L, angles)

    best = argmin(r -> r.shape_difference, result)
    @test isapprox(best.angle, 0.0; atol=1e-12)
    # comparing a shape to itself at zero rotation is an exact match
    @test isapprox(best.shape_difference, 0.0; atol=1e-10)
end

@testitem "shape_difference_rotation_boundary honours an injected metric" setup = [BoundaryRegSetup] begin
    called = Ref(0)
    counting_metric = function (a, b)
        called[] += 1
        return boundary_mse_aligned(a, b)
    end

    angles = [0.0, 0.1, 0.2]
    result = shape_difference_rotation_boundary(L, L, angles; metric=counting_metric)

    @test called[] == length(angles)          # metric actually used, once per angle
    # and it produced the injected metric's values, not the default's
    expected = [r.shape_difference for r in
                shape_difference_rotation_boundary(L, L, angles; metric=boundary_mse_aligned)]
    @test [r.shape_difference for r in result] == expected
end

@testitem "register_boundary recovers a known rotation" setup = [BoundaryRegSetup] begin
    θ = deg2rad(30.0)                    # on the default 5 degree grid
    target = rotate_boundary(L, θ)

    recovered = register_boundary(L, target)
    @test isapprox(recovered, θ; atol=1e-9)
end

@testitem "register_boundary sign convention matches the mask path" setup = [BoundaryRegSetup] begin
    # The mask-based shape_difference_rotation rotates the TARGET by -angle, so a
    # returned angle A means "target looks like reference rotated by A". Pin that
    # same direction here: if target is the reference rotated by +A, we get +A back
    # (not -A), which is what makes register_boundary a drop-in for register in
    # get_rotation_measurements.
    for deg in (10.0, 25.0, -15.0, -40.0)
        θ = deg2rad(deg)
        target = rotate_boundary(L, θ)
        @test isapprox(register_boundary(L, target), θ; atol=1e-9)
    end
end

@testitem "register_boundary returns an angle from test_angles" setup = [BoundaryRegSetup] begin
    target = rotate_boundary(L, deg2rad(17.3))   # deliberately off-grid
    angles = register_default_angles_rad

    recovered = register_boundary(L, target; test_angles=angles)
    @test any(a -> isapprox(a, recovered; atol=1e-12), angles)
end

@testitem "register_boundary works with a prior-restricted grid" setup = [BoundaryRegSetup] begin
    θ = deg2rad(12.0)
    target = rotate_boundary(L, θ)

    angles = prior_test_angles(θ; window=deg2rad(10.0))
    recovered = register_boundary(L, target; test_angles=angles)
    @test isapprox(recovered, θ; atol=1e-9)
end

@testitem "register_boundary is callable as a registration_function" setup = [BoundaryRegSetup] begin
    # get_rotation_measurements invokes registration_function(image1, image2; test_angles)
    θ = deg2rad(20.0)
    target = rotate_boundary(L, θ)

    f = (a, b; test_angles) -> register_boundary(a, b; test_angles)
    @test isapprox(f(L, target; test_angles=register_default_angles_rad), θ; atol=1e-9)
end
