@testitem "utils.jl pad utilities" begin
    using Images: Pad, Fill

    fill_value = 0
    simpleimg = collect(reshape(1:4, 2, 2))
    w, h = size(simpleimg)
    # 2×2 Matrix{Int64}:
    #  1  3
    #  2  4
    lo = (1, 2)
    hi = (3, 4)
    fil_type = Fill(fill_value, lo, hi)
    rep_type = Pad(:replicate, lo, hi)
    fil_paddedimg = IceFloeTracker.add_padding(simpleimg, fil_type)
    rep_paddedimg = IceFloeTracker.add_padding(simpleimg, rep_type)

    # test replicate padding at corners
    @test rep_paddedimg[1, 1] == 1 &&
        rep_paddedimg[1, end] == 3 &&
        rep_paddedimg[end, 1] == 2 &&
        rep_paddedimg[end, end] == 4

    # test filled image at the corners
    @test fil_paddedimg[1, 1] ==
        fil_paddedimg[1, end] ==
        fil_paddedimg[end, 1] ==
        fil_paddedimg[end, end] ==
        0

    # test sizing
    @test size(fil_paddedimg) == (h + lo[1] + hi[1], w + lo[2] + hi[2])

    # test padding removal
    @test IceFloeTracker.remove_padding(fil_paddedimg, fil_type) == simpleimg
    @test IceFloeTracker.remove_padding(rep_paddedimg, rep_type) == simpleimg
end
