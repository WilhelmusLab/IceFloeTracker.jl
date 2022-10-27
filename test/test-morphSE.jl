# Start with a pixel in the middle and dilate in one go to fill up the full image
@testset "MorphSE test" begin
@time begin
    n = rand(11:2:21) # choose random odd number
    mid = (n - 1) ÷ 2 + 1 # get median
    a = zeros(Int, n, n) # 
    a[mid, mid] = 1 # make 1 the pixel in the center
    se = IceFloeTracker.MorphSE.strel_box((n, n)) 
    @test IceFloeTracker.MorphSE.dilate(a, se) == ones(Int, n, n)
end
end