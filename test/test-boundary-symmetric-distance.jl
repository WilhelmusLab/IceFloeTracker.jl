# boundary_symmetric_distance: the mask path's centred symmetric difference
# (count_symdiff after align_centroids), evaluated on closed boundary curves by
# rasterising the regions they enclose. It exists to validate the boundary path
# against the mask path in the same units (pixels), not for tracker throughput.

@testsnippet SymDiffSetup begin
    using IceFloeTracker.Tracking:
        boundary_symmetric_distance,
        count_symdiff,
        align_centroids,
        rotate_boundary,
        register,
        register_boundary,
        imrotate_bin,
        bwtraceboundary,
        resample_boundary

    # closed, integer-cornered so the rasterised areas are exact
    square10 = [0.0 0.0; 10.0 0.0; 10.0 10.0; 0.0 10.0; 0.0 0.0]
    rect10x14 = [0.0 0.0; 10.0 0.0; 10.0 14.0; 0.0 14.0; 0.0 0.0]
    L = [0.0 0.0; 3.0 0.0; 3.0 1.0; 1.0 1.0; 1.0 3.0; 0.0 3.0; 0.0 0.0]

    # the same closed curve, traced from a different starting vertex
    function retrace(b, k)
        pts = b[1:(end - 1), :]
        rolled = circshift(pts, (-k, 0))
        return vcat(rolled, rolled[1:1, :])
    end

    # what add_boundary! does, without needing a DataFrame
    function traced(mask; reduc_factor=2)
        bd = bwtraceboundary(mask)
        bd1 = isa(bd, Vector{Vector{CartesianIndex}}) ? bd[1] : bd
        return resample_boundary(bd1, reduc_factor)
    end

    # an L in a frame wide enough that a 30° rotation stays inside it
    floe = falses(60, 60)
    floe[15:45, 15:30] .= true
    floe[15:25, 30:45] .= true
end

@testitem "boundary_symmetric_distance: identity is zero" setup = [SymDiffSetup] begin
    @test boundary_symmetric_distance(L, L) == 0
    @test boundary_symmetric_distance(square10, square10) == 0
end

@testitem "boundary_symmetric_distance: units are pixels, exact on lattice-aligned shapes" setup = [SymDiffSetup] begin
    # 10x10 against 10x14, both centred: the extra 4 rows are 40 pixels, all of
    # them in the symmetric difference. Same integer count_symdiff would return.
    @test boundary_symmetric_distance(square10, rect10x14) == 40
end

@testitem "boundary_symmetric_distance: symmetric" setup = [SymDiffSetup] begin
    @test boundary_symmetric_distance(L, rect10x14) == boundary_symmetric_distance(rect10x14, L)
end

@testitem "boundary_symmetric_distance: invariant to trace starting point" setup = [SymDiffSetup] begin
    for k in 1:(size(L, 1) - 2)
        @test boundary_symmetric_distance(L, retrace(L, k)) == 0
    end
end

@testitem "boundary_symmetric_distance: invariant to translation" setup = [SymDiffSetup] begin
    # centred internally, like align_centroids on the mask path
    @test boundary_symmetric_distance(L, L .+ [17.0 -4.0]) == 0
    @test boundary_symmetric_distance(square10, square10 .+ [0.5 0.25]) == 0
end

@testitem "boundary_symmetric_distance: sees rotation" setup = [SymDiffSetup] begin
    # rotation is what register_boundary searches over, so it must not cancel
    @test boundary_symmetric_distance(L, rotate_boundary(L, deg2rad(45.0))) > 0
end

@testitem "boundary_symmetric_distance agrees with the mask path on a traced floe" setup = [SymDiffSetup] begin
    # The reason this metric exists. The boundary is a resampled B-spline through
    # the pixel ring, so its rasterisation is not the mask bit for bit; the
    # agreement asked of it is within a one-pixel band along the perimeter.
    θ = deg2rad(30.0)
    m1, m2 = floe, imrotate_bin(floe, θ)
    b1, b2 = traced(m1), traced(m2)

    mask_value = count_symdiff(align_centroids(m1, m2)...)
    bnd_value = boundary_symmetric_distance(b1, b2)

    ring = bwtraceboundary(m1)
    perimeter_px = length(isa(ring, Vector{Vector{CartesianIndex}}) ? ring[1] : ring)
    @test abs(bnd_value - mask_value) <= perimeter_px
    @test bnd_value > 0                         # a 30° rotation is not a no-op
end

@testitem "register_boundary with boundary_symmetric_distance returns register's angle" setup = [SymDiffSetup] begin
    # Exact comparability, as the reviewer put it: the same objective the mask
    # path minimises, evaluated on the curve. Both paths must then land on the
    # same grid angle for the same physical rotation, sign included.
    θ = deg2rad(30.0)                                  # on the default 5° grid
    m1, m2 = floe, imrotate_bin(floe, θ)
    b1, b2 = traced(m1), traced(m2)

    from_mask = register(m1, m2)
    from_boundary = register_boundary(b1, b2; metric=boundary_symmetric_distance)

    @test isapprox(from_mask, θ; atol=1e-9)
    @test isapprox(from_boundary, from_mask; atol=1e-9)
end
