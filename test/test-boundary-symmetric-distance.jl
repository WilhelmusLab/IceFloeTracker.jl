# boundary_symmetric_distance: the mask path's centered symmetric difference
# (count_symdiff after align_centroids), evaluated on closed boundary curves by
# rasterizing the regions they enclose. It exists to validate the boundary path
# against the mask path in the same units (pixels), not for tracker throughput.

@testsnippet SymDiffSetup begin
    using IceFloeTracker.Tracking:
        boundary_symmetric_distance,
        _polygon_centroid,
        count_symdiff,
        align_centroids,
        rotate_boundary,
        register,
        register_boundary,
        imrotate_bin,
        bwtraceboundary,
        resample_boundary

    # closed, integer-cornered so the rasterized areas are exact
    square10 = [0.0 0.0; 10.0 0.0; 10.0 10.0; 0.0 10.0; 0.0 0.0]
    rect10x14 = [0.0 0.0; 10.0 0.0; 10.0 14.0; 0.0 14.0; 0.0 0.0]
    L = [0.0 0.0; 3.0 0.0; 3.0 1.0; 1.0 1.0; 1.0 3.0; 0.0 3.0; 0.0 0.0]
    #     L  (6 vertices)
    # 3 o───────o
    #   │       │
    #   │       │
    #   │       │
    # 2 │       │
    #   │       │
    #   │       │
    #   │       │  +
    # 1 │       o×──────────────o
    #   │                       │
    #   │                       │
    #   │                       │
    # 0 o───────────────────────o
    #   0       1       2       3

    # The L with collinear vertices added on two of its edges: same polygon, same
    # centroid, while the vertex mean moves.
    # dense  (12 vertices)
    # 3 o───o───o
    #   │       │
    #   │       │
    #   │       │
    # 2 │       │
    #   │       │
    #   │       │
    #   │       │
    # 1 │       o×─+────────────o
    #   │                       │
    #   │                       │
    #   │                       │
    # 0 o───o───o───o───o───o───o
    #   0       1       2       3
    extra_vertices = [0.5 0.0; 1.0 0.0; 1.5 0.0; 2.0 0.0; 2.5 0.0] # along bottom
    Ldense = vcat(
        L[1:1, :], extra_vertices, L[2:5, :], [0.5 3.0], L[6:end, :] # one more on the top
    )

    # the same closed curve, traced from a different starting vertex
    function retrace(b, k)
        pts = b[1:(end-1), :]
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
    # 10x10 against 10x14, both centered: the extra 4 rows are 40 pixels, all of
    # them in the symmetric difference. Same integer count_symdiff would return.
    @test boundary_symmetric_distance(square10, rect10x14) == 40
end

@testitem "boundary_symmetric_distance: edges on pixel centers, either orientation" setup = [SymDiffSetup] begin
    using IceFloeTracker.Tracking: _rasterize
    # odd sides put the centered edges on pixel centers; each pixel must be counted once,
    # whichever axis the edge lies along. Once centered, the pixel centers fall on the
    # integer points below (■ filled, □ not filled; the corners are the matrix rows):
    #
    # rect3x5               rect5x3
    #  y                     y
    #  5  □───□───□───□      3  □───□───□───□───□───□
    #     │           │         │                   │
    #  4  ■   ■   ■   □      2  ■   ■   ■   ■   ■   □
    #     │           │         │                   │
    #  3  ■   ■   ■   □      1  ■   ■   ■   ■   ■   □
    #     │           │         │                   │
    #  2  ■   ■   ■   □      0  ■───■───■───■───■───□
    #     │           │         0   1   2   3   4   5  x
    #  1  ■   ■   ■   □
    #     │           │
    #  0  ■───■───■───□
    #     0   1   2   3  x
    #
    # Left and bottom edges are filled, right and top are not
    rect3x5 = [0.0 0.0; 3.0 0.0; 3.0 5.0; 0.0 5.0; 0.0 0.0]
    rect5x3 = rect3x5[:, [2, 1]]
    for r in (rect3x5, rect5x3)
        c = r .- _polygon_centroid(r)'
        @test count(_rasterize(c, (-3, -3), (3, 3))) == 15
    end
    # centered overlap is 3x3, so 15 + 15 - 2 * 9
    @test boundary_symmetric_distance(rect3x5, rect5x3) == 12
end

@testitem "boundary_symmetric_distance: symmetric" setup = [SymDiffSetup] begin
    @test boundary_symmetric_distance(L, rect10x14) == boundary_symmetric_distance(rect10x14, L)
end

@testitem "boundary_symmetric_distance: invariant to trace starting point" setup = [SymDiffSetup] begin
    for k in 1:(size(L, 1)-2)
        @test boundary_symmetric_distance(L, retrace(L, k)) == 0
    end
end

@testitem "boundary_symmetric_distance: invariant to translation" setup = [SymDiffSetup] begin
    # centered internally, like align_centroids on the mask path
    @test boundary_symmetric_distance(L, L .+ [17.0 -4.0]) == 0
    @test boundary_symmetric_distance(square10, square10 .+ [0.5 0.25]) == 0
end

@testitem "_polygon_centroid is the area centroid, not the vertex mean" setup = [SymDiffSetup] begin
    # L = 3x1 bar (centroid (1.5, 0.5)) + 1x2 bar (centroid (0.5, 2.0)), area 5:
    # (1.1, 1.1). The mean of its six vertices is (4/3, 4/3).
    @test all(isapprox.(_polygon_centroid(L), [1.1, 1.1]; atol=1e-12))
    @test all(isapprox.(_polygon_centroid(reverse(L; dims=1)), [1.1, 1.1]; atol=1e-12))
    @test all(isapprox.(_polygon_centroid(L[1:(end-1), :]), [1.1, 1.1]; atol=1e-12))
end

@testitem "_polygon_centroid is invariant to extra vertices along an edge" setup = [SymDiffSetup] begin
    using Statistics: mean
    # Check vertex mean moves
    @test !isapprox(vec(mean(Ldense[1:(end-1), :]; dims=1)), [1.1, 1.1]; atol=1e-3)

    # Check that the area centroid remains the same
    @test all(isapprox.(_polygon_centroid(Ldense), _polygon_centroid(L); atol=1e-12))
end

@testitem "boundary_symmetric_distance: invariant to vertex density" setup = [SymDiffSetup] begin
    # Ldense scaled so the vertex-mean offset spans several pixels; at unit scale a
    # sub-pixel offset can leave the rasterized count unchanged.
    @test boundary_symmetric_distance(10 .* L, 10 .* Ldense) == 0
end

@testitem "boundary_symmetric_distance: closure only up to round-off" setup = [SymDiffSetup] begin
    # Resampled boundaries close only to ~1e-13, so the closing row must not be
    # counted as an extra vertex.
    L10 = 10 .* L
    for k in 0:(size(L10, 1)-2)
        r = retrace(L10, k)
        r[end, :] .+= 1e-13
        @test boundary_symmetric_distance(L10, r) == 0
    end
end

@testitem "boundary_symmetric_distance: sees rotation" setup = [SymDiffSetup] begin
    # rotation is what register_boundary searches over, so it must not cancel
    @test boundary_symmetric_distance(L, rotate_boundary(L, deg2rad(45.0))) > 0
end

@testitem "boundary_symmetric_distance agrees with the mask path on a traced floe" setup = [SymDiffSetup] begin
    # The reason this metric exists. The boundary is a resampled B-spline through
    # the pixel ring, so its rasterization is not the mask bit for bit; the
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
    # path minimizes, evaluated on the curve. Both paths must then land on the
    # same grid angle for the same physical rotation, sign included.
    θ = deg2rad(30.0)                                  # on the default 5° grid
    m1, m2 = floe, imrotate_bin(floe, θ)
    b1, b2 = traced(m1), traced(m2)

    from_mask = register(m1, m2)
    from_boundary = register_boundary(b1, b2; metric=boundary_symmetric_distance)

    @test isapprox(from_mask, θ; atol=1e-9)
    @test isapprox(from_boundary, from_mask; atol=1e-9)
end
