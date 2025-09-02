
@testitem "imcomplement" begin
    using IceFloeTracker: imcomplement
    using Images: Gray

    img = rand(0:255, 10, 10)
    @test imcomplement(img) == 255 .- img

    img = Gray.(img ./ 255)
    @test imcomplement(img) == 1 .- img
end
