@testitem "BenkridCrookes perimeter" begin
    using IceFloeTracker: BenkridCrookes

    bc4 = BenkridCrookes(; connectivity=4)
    bc8 = BenkridCrookes(; connectivity=8)

    # The two cases in the `BenkridCrookes` docstring.
    A = [0 1 1; 1 1 1; 1 1 1]
    @test bc4(A) == 7.414213562373095
    @test bc8(A) == 7.0

    # A lone pixel has no neighborhood to trace around.
    @test bc4(trues(1, 1)) == 0.0

    # The estimate follows the boundary pixel centers, so a solid k x k square
    # measures 4(k - 1) rather than 4k.
    for k in 2:8
        @test bc4(trues(k, k)) == 4.0 * (k - 1)
    end

    # A shape and its transpose, and a shape and its mirror, enclose the same
    # perimeter.
    shapes = [
        Bool[1 1 0; 1 1 1; 0 1 1],
        Bool[1 0 0 0; 1 1 1 0; 0 1 1 1],
        Bool[0 1 0; 1 1 1; 0 1 0],
        trues(3, 7),
    ]
    for s in shapes
        @test bc4(s) == bc4(permutedims(s))
        @test bc4(s) == bc4(reverse(s; dims=1))
        @test bc4(s) == bc4(reverse(s; dims=2))
        @test bc8(s) == bc8(permutedims(s))
    end

    # Out of frame counts as background, so surrounding a mask with background
    # cannot change its perimeter. Masks arrive cropped to their bounding box,
    # which makes the foreground-on-the-frame-edge case the common one.
    for s in shapes
        n, m = size(s)
        padded = falses(n + 4, m + 4)
        padded[3:(n+2), 3:(m+2)] .= s
        @test bc4(s) == bc4(padded)
        @test bc8(s) == bc8(padded)
    end

    # Masks reach this function as a BitMatrix, the docstring passes a Matrix{Int}:
    # anything where zero means background has to agree.
    @test bc4(Bool.(A)) == bc4(A)
    @test bc4(BitMatrix(Bool.(A))) == bc4(A)
end

@testitem "BenkridCrookes boundary image" begin
    using IceFloeTracker.Segmentation: _boundary_image

    # The result is zero-padded by one pixel, so callers can read
    # every pixel's neighborhood without a bounds test.
    e = _boundary_image(trues(3, 4), true)
    @test size(e) == (5, 6)
    @test !any(e[1, :]) && !any(e[end, :]) && !any(e[:, 1]) && !any(e[:, end])

    # A solid square is boundary everywhere except its interior.
    e = _boundary_image(trues(3, 3), true)
    @test e[2:4, 2:4] == Bool[1 1 1; 1 0 1; 1 1 1]

    # A lone pixel is all boundary; empty input is all background.
    @test _boundary_image(trues(1, 1), true)[2, 2]
    @test !any(_boundary_image(falses(4, 4), true))

    # Connectivity decides whether diagonal neighbors count. The center of a
    # plus has all four edge neighbors but no diagonal ones, so it is interior
    # under 4-connectivity and boundary under 8.
    plus = Bool[0 1 0; 1 1 1; 0 1 0]
    @test !_boundary_image(plus, true)[3, 3]
    @test _boundary_image(plus, false)[3, 3]

    # Zero means background whatever the element type.
    @test _boundary_image([0 1; 1 1], true) == _boundary_image(Bool[0 1; 1 1], true)
end

@testitem "BenkridCrookes neighborhood codes" begin
    using IceFloeTracker.Segmentation: _boundary_image, _boundary_type_counts

    # A code is 1 + 2*(edge neighbors on the boundary) + 10*(diagonal ones), so
    # it runs to 1 + 2*4 + 10*4 = 49.
    counts = _boundary_type_counts(_boundary_image(trues(2, 2), true))
    @test length(counts) == 49

    # Each pixel of a 2x2 square sees two edge neighbors and one diagonal:
    # 1 + 2*2 + 10*1 = 15.
    @test counts[15] == 4
    @test sum(counts) == 4

    # A lone pixel has an empty neighborhood: 1 + 0 + 0.
    @test _boundary_type_counts(_boundary_image(trues(1, 1), true))[1] == 1

    # Nothing to classify.
    @test all(iszero, _boundary_type_counts(_boundary_image(falses(3, 3), true)))

    # Every boundary pixel is counted exactly once.
    shapes = (Bool[1 1 0; 1 1 1; 0 1 1], trues(4, 6), Bool[1 0 1; 0 1 0; 1 0 1])
    for s in shapes
        e = _boundary_image(s, true)
        @test sum(_boundary_type_counts(e)) == count(e)
    end
end
