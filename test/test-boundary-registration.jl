# Phase 4: boundary-curve analogues of shape_difference_rotation / register.
#
# Fixtures here are CLOSED (first point repeated as last), matching what
# add_boundary! actually produces. They are also ASYMMETRIC under rotation --
# a square would make the recovered angle ambiguous at 90 degree multiples.

@testsnippet BoundaryRegSetup begin
    using IceFloeTracker.Tracking:
        shape_difference_rotation_boundary,
        register_boundary,
        boundary_shape_difference,
        BoundaryRegistration,
        rotate_boundary,
        boundary_hausdorff,
        boundary_modified_hausdorff,
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

    # The boundary that tracing imrotate_bin(mask, θ) would give, up to pixelisation.
    # Boundary coordinates are image indices (row axis down), so imrotate's positive
    # sense is rotate_boundary's negative. Targets are built with this so the angle
    # register_boundary must return is the same one register returns for the masks.
    rotate_as_mask(b, θ) = rotate_boundary(b, -θ)
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
        return boundary_hausdorff(a, b)
    end

    angles = [0.0, 0.1, 0.2]
    result = shape_difference_rotation_boundary(L, L, angles; metric=counting_metric)

    @test called[] == length(angles)          # metric actually used, once per angle
    # and it produced the injected metric's values, not the default's
    expected = [r.shape_difference for r in
                shape_difference_rotation_boundary(L, L, angles; metric=boundary_hausdorff)]
    @test [r.shape_difference for r in result] == expected
end

@testitem "register_boundary recovers a known rotation" setup = [BoundaryRegSetup] begin
    θ = deg2rad(30.0)                    # on the default 5 degree grid
    target = rotate_as_mask(L, θ)

    recovered = register_boundary(L, target)
    @test isapprox(recovered, θ; atol=1e-9)
end

@testitem "register_boundary sign convention matches the mask path" setup = [BoundaryRegSetup] begin
    # A returned angle A means "target looks like reference rotated by A" in the sense
    # of imrotate_bin_clockwise_radians, the convention register uses. Pinning this is
    # what makes register_boundary a drop-in for register in get_rotation_measurements;
    # the traced-boundary test below checks it against register itself.
    for deg in (10.0, 25.0, -15.0, -40.0)
        θ = deg2rad(deg)
        target = rotate_as_mask(L, θ)
        @test isapprox(register_boundary(L, target), θ; atol=1e-9)
    end
end

@testitem "register_boundary returns an angle from test_angles" setup = [BoundaryRegSetup] begin
    target = rotate_as_mask(L, deg2rad(17.3))   # deliberately off-grid
    angles = register_default_angles_rad

    recovered = register_boundary(L, target; test_angles=angles)
    @test any(a -> isapprox(a, recovered; atol=1e-12), angles)
end

@testitem "register_boundary works with a prior-restricted grid" setup = [BoundaryRegSetup] begin
    θ = deg2rad(12.0)
    target = rotate_as_mask(L, θ)

    angles = prior_test_angles(θ; window=deg2rad(10.0))
    recovered = register_boundary(L, target; test_angles=angles)
    @test isapprox(recovered, θ; atol=1e-9)
end

@testitem "register_boundary recovers rotation from independently traced boundaries" setup = [BoundaryRegSetup] begin
    # The production case. The two boundaries are traced from separately produced masks,
    # so their point sequences start at unrelated places and carry pixel noise. Every
    # other test here builds its target with rotate_boundary, which preserves point
    # order and so cannot detect a metric that depends on it -- an index-wise metric
    # passed all of those and recovered the angle on 2% of real floes.
    using IceFloeTracker.Tracking: _traced_boundary, imrotate_bin

    # an L: no rotational symmetry, and large enough that a 30° rotation survives
    # pixelisation without clipping inside the 60x60 frame
    mask = falses(60, 60)
    mask[15:45, 15:30] .= true
    mask[15:25, 30:45] .= true

    θ = deg2rad(30.0)                                    # on the default 5° grid
    reference = _traced_boundary(mask)
    target = _traced_boundary(imrotate_bin(mask, θ))    # traced anew, not rotated analytically

    using IceFloeTracker.Tracking: register
    recovered = register_boundary(reference, target)
    @test isapprox(recovered, θ; atol=1e-9)
    @test isapprox(recovered, register(mask, imrotate_bin(mask, θ)); atol=1e-9)
end

@testitem "register_boundary is callable as a registration_function" setup = [BoundaryRegSetup] begin
    # get_rotation_measurements invokes registration_function(image1, image2; test_angles)
    θ = deg2rad(20.0)
    target = rotate_as_mask(L, θ)

    f = (a, b; test_angles) -> register_boundary(a, b; test_angles)
    @test isapprox(f(L, target; test_angles=register_default_angles_rad), θ; atol=1e-9)
end

# ---------------------------------------------------------------------------
# BoundaryRegistration: a configured, callable wrapper around register_boundary
# ---------------------------------------------------------------------------

@testitem "BoundaryRegistration defaults" setup = [BoundaryRegSetup] begin
    reg = BoundaryRegistration()
    @test reg.metric === boundary_modified_hausdorff
    @test reg.boundary_column === :boundary

    # Must subtype Function or rotation.jl's `registration_function::Function`
    # keyword rejects it with a TypeError at the call site.
    @test reg isa Function
end

@testitem "BoundaryRegistration matrix method matches register_boundary" setup = [BoundaryRegSetup] begin
    θ = deg2rad(30.0)
    target = rotate_as_mask(L, θ)
    reg = BoundaryRegistration()

    @test reg(L, target) == register_boundary(L, target)   # same grid, exact
    @test isapprox(reg(L, target), θ; atol=1e-9)
    @test isfinite(reg(L, target))
end

@testitem "BoundaryRegistration honours test_angles" setup = [BoundaryRegSetup] begin
    θ = deg2rad(12.0)
    target = rotate_as_mask(L, θ)
    angles = prior_test_angles(θ; window=deg2rad(10.0))

    @test isapprox(BoundaryRegistration()(L, target; test_angles=angles), θ; atol=1e-9)
end

@testitem "BoundaryRegistration uses the injected metric" setup = [BoundaryRegSetup] begin
    called = Ref(0)
    counting = function (a, b)
        called[] += 1
        return boundary_hausdorff(a, b)
    end

    θ = deg2rad(20.0)
    target = rotate_as_mask(L, θ)
    angles = [deg2rad(d) for d in 15:25]

    got = BoundaryRegistration(; metric=counting)(L, target; test_angles=angles)

    @test called[] == length(angles)   # metric invoked once per angle
    @test got == register_boundary(L, target; test_angles=angles, metric=boundary_hausdorff)
end

@testitem "BoundaryRegistration all metrics recover the angle" setup = [BoundaryRegSetup] begin
    θ = deg2rad(25.0)
    target = rotate_as_mask(L, θ)
    for m in (boundary_modified_hausdorff, boundary_hausdorff)
        @test isapprox(BoundaryRegistration(; metric=m)(L, target), θ; atol=1e-9)
    end
end

@testitem "BoundaryRegistration DataFrameRow method" setup = [BoundaryRegSetup] begin
    using DataFrames

    θ = deg2rad(25.0)
    df = DataFrame(; boundary=[L, rotate_as_mask(L, θ)])
    reg = BoundaryRegistration()

    @test isapprox(reg(df[1, :], df[2, :]), θ; atol=1e-9)
    @test reg(df[1, :], df[2, :]) == reg(L, rotate_as_mask(L, θ))
    # kwargs must splat through the row method
    @test isapprox(
        reg(df[1, :], df[2, :]; test_angles=register_default_angles_rad), θ; atol=1e-9
    )

    # a DataFrameRow must not be an AbstractMatrix, or it would dispatch to the
    # widening matrix method instead of the row method
    @test !(typeof(df[1, :]) <: AbstractMatrix)
end

@testitem "BoundaryRegistration configurable boundary_column" setup = [BoundaryRegSetup] begin
    using DataFrames

    θ = deg2rad(15.0)
    df = DataFrame(; bd=[L, rotate_as_mask(L, θ)])
    reg = BoundaryRegistration(; boundary_column=:bd)

    @test isapprox(reg(df[1, :], df[2, :]), θ; atol=1e-9)
end

@testitem "BoundaryRegistration widens non-Float64 boundaries" setup = [BoundaryRegSetup] begin
    # the :boundary column is untyped, so an Int matrix can reach the functor
    square = [0 0; 3 0; 3 1; 1 1; 1 3; 0 3; 0 0]
    @test eltype(square) == Int
    @test isapprox(BoundaryRegistration()(square, square), 0.0; atol=1e-9)
end

# ---------------------------------------------------------------------------
# boundary_shape_difference: orientation-aligned, non-searching comparison
# (boundary analogue of shape_difference, which the candidate filter needs)
# ---------------------------------------------------------------------------

@testitem "boundary_shape_difference identity and symmetry" setup = [BoundaryRegSetup] begin
    @test isapprox(boundary_shape_difference(L, 0.0, L, 0.0), 0.0; atol=1e-10)
    @test isapprox(
        boundary_shape_difference(L, 0.3, L, 0.3), 0.0; atol=1e-10
    )  # same rotation applied to both

    other = rotate_as_mask(L, deg2rad(40.0))
    @test isapprox(
        boundary_shape_difference(L, 0.1, other, 0.2),
        boundary_shape_difference(other, 0.2, L, 0.1);
        atol=1e-10,
    )
end

@testitem "boundary_shape_difference aligns by orientation" setup = [BoundaryRegSetup] begin
    # a shape and its rotation, each labelled with its own orientation, should
    # align to near-zero difference once each is rotated by that orientation
    θ = deg2rad(35.0)
    rotated = rotate_as_mask(L, θ)
    @test boundary_shape_difference(L, 0.0, rotated, -θ) <
        boundary_shape_difference(L, 0.0, rotated, 0.0)
end

@testitem "boundary_shape_difference honours the metric kwarg" setup = [BoundaryRegSetup] begin
    other = rotate_as_mask(L, deg2rad(40.0))
    a = boundary_shape_difference(L, 0.0, other, 0.0; metric=boundary_hausdorff)
    b = boundary_shape_difference(L, 0.0, other, 0.0; metric=boundary_modified_hausdorff)
    @test a != b            # max vs mean nearest-neighbour distance differ unless all are equal
    @test isfinite(a) && isfinite(b)
end

@testitem "boundary_shape_difference DataFrameRow method" setup = [BoundaryRegSetup] begin
    using DataFrames

    df = DataFrame(; boundary=[L, L], orientation=[0.0, 0.0])
    @test isapprox(boundary_shape_difference(df[1, :], df[2, :]), 0.0; atol=1e-10)
end

# ---------------------------------------------------------------------------
# End-to-end: BoundaryRegistration as a registration_function.
# These are the only tests that exercise rotation.jl's injection path, and the
# prior variant is the only one that exercises its `test_angles` closure.
# ---------------------------------------------------------------------------

@testitem "BoundaryRegistration drives get_rotation_measurements" setup = [BoundaryRegSetup] begin
    using DataFrames, Dates
    using IceFloeTracker: get_rotation_measurements

    θ = deg2rad(30.0)                      # on the default 5 degree grid
    t0 = DateTime(2020, 1, 1, 0, 0, 0)

    # row1 is the older observation; prior = orientation1 - orientation2, so
    # orientation2 = -θ makes the prior +θ (matching the mask-based convention)
    df = DataFrame(;
        id=[1, 1],
        boundary=[L, rotate_as_mask(L, θ)],
        orientation=[0.0, -θ],
        time=[t0, t0 + Hour(6)],
    )

    result = get_rotation_measurements(
        df;
        id_column=:id,
        image_column=:boundary,
        time_column=:time,
        registration_function=BoundaryRegistration(),
    )

    @test nrow(result) == 1
    @test result.theta_rad[1] isa Float64
    @test isapprox(result.theta_rad[1], θ; atol=1e-9)
    @test result.dt_sec[1] == 6 * 3600
    @test isapprox(result.omega_rad_per_sec[1], θ / (6 * 3600); atol=1e-12)
    @test result.omega_rad_per_sec[1] > 0
end

@testitem "BoundaryRegistration end-to-end with a negative rotation" setup = [BoundaryRegSetup] begin
    using DataFrames, Dates
    using IceFloeTracker: get_rotation_measurements

    θ = deg2rad(-25.0)
    t0 = DateTime(2020, 1, 1, 0, 0, 0)
    df = DataFrame(;
        id=[1, 1],
        boundary=[L, rotate_as_mask(L, θ)],
        orientation=[0.0, -θ],
        time=[t0, t0 + Hour(6)],
    )

    result = get_rotation_measurements(
        df;
        id_column=:id,
        image_column=:boundary,
        time_column=:time,
        registration_function=BoundaryRegistration(),
    )

    @test isapprox(result.theta_rad[1], θ; atol=1e-9)
    @test result.omega_rad_per_sec[1] < 0        # sign must survive the pipeline
end

@testitem "BoundaryRegistration end-to-end honours the orientation prior" setup = [BoundaryRegSetup] begin
    using DataFrames, Dates
    using IceFloeTracker: get_rotation_measurements

    θ = deg2rad(30.0)
    t0 = DateTime(2020, 1, 1, 0, 0, 0)
    df = DataFrame(;
        id=[1, 1],
        boundary=[L, rotate_as_mask(L, θ)],
        orientation=[0.0, -θ],
        time=[t0, t0 + Hour(6)],
    )
    kw = (;
        id_column=:id,
        image_column=:boundary,
        time_column=:time,
        registration_function=BoundaryRegistration(),
    )

    # guard: if the prior grid does not contain the true angle, a failure below
    # would be about prior construction, not about BoundaryRegistration
    @test any(a -> isapprox(a, θ; atol=1e-9), prior_test_angles(θ; window=deg2rad(15.0)))

    full = get_rotation_measurements(df; kw...)
    with_prior = get_rotation_measurements(
        df; kw..., orientation_column=:orientation, orientation_window=deg2rad(15.0)
    )

    # Unlike the mask path -- whose equivalents in test-rotation.jl are marked
    # broken because pixelated imrotate on a fine grid cannot reproduce the
    # coarse-grid optimum -- boundary rotation is exact, so these agree.
    @test isapprox(with_prior.theta_rad[1], full.theta_rad[1]; atol=1e-8)
    @test isapprox(with_prior.theta_rad[1], θ; atol=1e-9)
end

@testitem "BoundaryRegistration survives an empty DataFrame" setup = [BoundaryRegSetup] begin
    using DataFrames, Dates
    using IceFloeTracker: get_rotation_measurements

    df = DataFrame(;
        id=Int[], boundary=Matrix{Float64}[], orientation=Float64[], time=DateTime[]
    )
    result = get_rotation_measurements(
        df;
        id_column=:id,
        image_column=:boundary,
        time_column=:time,
        registration_function=BoundaryRegistration(),
    )
    @test nrow(result) == 0
end
