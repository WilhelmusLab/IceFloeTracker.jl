@testitem "imadjust" begin
    using Random
    Random.seed!(123)
    img = rand(0:255, 100, 100)
    @test sum(IceFloeTracker.imadjust(img)) == 1291155
end
