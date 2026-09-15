# Set-based boundary metrics: Hausdorff and modified Hausdorff.
#
# Fixtures are CLOSED (first row repeated as last), matching add_boundary! output, and
# ASYMMETRIC so that rotations are distinguishable.

@testsnippet SetMetricSetup begin
    using IceFloeTracker.Tracking:
        boundary_hausdorff, boundary_modified_hausdorff, rotate_boundary, center_boundary

    L = [
        0.0 0.0
        3.0 0.0
        3.0 1.0
        1.0 1.0
        1.0 3.0
        0.0 3.0
        0.0 0.0
    ]
    BAR = [
        0.0 0.0
        3.0 0.0
        3.0 0.4
        0.0 0.4
        0.0 0.0
    ]

    # the same closed curve, traced from a different starting point
    function retrace(b, k)
        open_pts = b[1:(end - 1), :]
        rolled = circshift(open_pts, (-k, 0))
        return vcat(rolled, rolled[1:1, :])
    end

    METRICS = (boundary_hausdorff, boundary_modified_hausdorff)
end

@testitem "set metrics: identity is zero" setup = [SetMetricSetup] begin
    for m in METRICS
        @test isapprox(m(L, L), 0.0; atol=1e-12)
    end
end

@testitem "set metrics: symmetric" setup = [SetMetricSetup] begin
    for m in METRICS
        @test isapprox(m(L, BAR), m(BAR, L); atol=1e-12)
    end
end

@testitem "set metrics: invariant to trace starting point" setup = [SetMetricSetup] begin
    # This is the property the previous metrics lacked. An index-wise comparison of L
    # against retrace(L, 2) scored 9.2 -- larger than L against a genuinely different
    # shape. A set metric must score the identical shape as identical.
    for m in METRICS, k in 1:(size(L, 1) - 2)
        @test isapprox(m(L, retrace(L, k)), 0.0; atol=1e-12)
    end
end

@testitem "set metrics: discriminate different shapes" setup = [SetMetricSetup] begin
    for m in METRICS
        @test m(L, BAR) > 0.1
    end
end

@testitem "set metrics: invariant to translation" setup = [SetMetricSetup] begin
    # both are centred internally, so a shifted copy is the same shape
    shifted = L .+ [17.0 -4.0]
    for m in METRICS
        @test isapprox(m(L, shifted), 0.0; atol=1e-12)
    end
end

@testitem "set metrics: not invariant to rotation" setup = [SetMetricSetup] begin
    # deliberately: rotation is what register_boundary searches over, so the metric
    # must see it
    for m in METRICS
        @test m(L, rotate_boundary(L, deg2rad(40.0))) > 0.1
    end
end

@testitem "modified Hausdorff is bounded by Hausdorff" setup = [SetMetricSetup] begin
    # a mean of nearest-neighbour distances cannot exceed their maximum
    for (a, b) in ((L, BAR), (L, rotate_boundary(L, 0.7)), (BAR, rotate_boundary(BAR, 1.1)))
        @test boundary_modified_hausdorff(a, b) <= boundary_hausdorff(a, b) + 1e-12
    end
end

@testitem "modified Hausdorff is robust to a single outlier" setup = [SetMetricSetup] begin
    # Dubuisson & Jain's motivation: one stray point (a segmentation artefact) should
    # not dominate the score. Perturb one vertex far out and compare the two responses.
    outlier = copy(L)
    outlier[3, :] .= [30.0, 30.0]     # move one interior vertex ~40 units away

    h = boundary_hausdorff(L, outlier)
    mhd = boundary_modified_hausdorff(L, outlier)

    @test h > 20.0                     # Hausdorff is set entirely by the outlier
    @test mhd < h / 3                  # the mean dilutes it by roughly 1/n
end

@testitem "set metrics: units are pixels, not dimensionless" setup = [SetMetricSetup] begin
    # scaling both shapes by s scales the score by s -- a length, unlike the
    # dimensionless ratio the previous default produced
    s = 4.0
    for m in METRICS
        @test isapprox(m(s .* L, s .* BAR), s * m(L, BAR); rtol=1e-9)
    end
end
