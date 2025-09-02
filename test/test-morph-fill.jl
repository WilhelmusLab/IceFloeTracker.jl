using IceFloeTracker: morph_fill

@testitem "Morph Fill" begin
    using IceFloeTracker: morph_fill

    img = Bool[
        0 0 0 0 0
        0 1 1 1 0
        0 1 0 1 0
        0 1 1 1 0
        0 0 0 0 0
    ]

    r = rand(1:100)

    imgbig = repeat(img, r, r)

    filled = morph_fill(imgbig)
    @test sum(filled) - sum(imgbig) == r^2
end
