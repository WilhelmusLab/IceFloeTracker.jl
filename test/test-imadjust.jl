@testset "imadjust" begin
    Random.seed!(123)
    img = rand(0:255, 100, 100)
    sum(IceFloeTracker.imadjust(img)) == 1291155
end