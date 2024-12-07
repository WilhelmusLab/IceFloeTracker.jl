using IceFloeTracker: imcomplement

@ntestset "$(@__FILE__)" begin
    img = rand(0:255, 10, 10)
    @test imcomplement(img) == 255 .- img

    img = Gray.(img ./ 255)
    @test imcomplement(img) == 1 .- img
end
