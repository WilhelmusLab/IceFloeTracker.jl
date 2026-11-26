@testitem "imhmin" begin
    import Images: erode, distance_transform, feature_transform

    # Reference: geodesic erosion of `img .+ h` above `img`, iterated until stable.
    function iterated_hmin(img, h)
        out = img .+ h
        while true
            next = max.(img, erode(out))
            next == out && return out
            out = next
        end
    end

    function disc(n, radius)
        center = (n + 1) / 2
        return [hypot(i - center, j - center) <= radius for i in 1:n, j in 1:n]
    end

    negative_distances(mask) = -distance_transform(feature_transform(.!mask))

    @testset "suppresses minima shallower than h" begin
        # Example from the MATLAB `imhmin` documentation
        img = fill(10.0, 10, 10)
        img[2:4, 2:4] .= 7 # depth 3: suppressed
        img[6:8, 6:8] .= 2 # depth 8: raised by h
        expected = fill(10.0, 10, 10)
        expected[6:8, 6:8] .= 6
        @test imhmin(img, 4) == expected
    end

    @testset "matches iterated geodesic erosion" begin
        floes = disc(41, 15) .|| circshift(disc(41, 8), (12, -14))
        images = (rand(30, 40), 10 .* rand(25, 25), negative_distances(floes))
        for img in images, h in (0.1, 0.5, 2)
            @test imhmin(img, h) == iterated_hmin(img, h)
        end
    end

    @testset "faster than iterated geodesic erosion" begin
        distances = negative_distances(disc(301, 140))
        seconds(f) = (f(); minimum(@elapsed(f()) for _ in 1:3))
        fast() = imhmin(distances, 2)
        reference() = iterated_hmin(distances, 2)
        # The reference makes one full-image pass per pixel of floe radius. Reconstruction
        # in a single pass measured ~20x faster on this floe, ~26x with coverage and bounds
        # checks on; 10x leaves margin for timing noise on CI.
        @test 10 * seconds(fast) < seconds(reference)
    end
end
