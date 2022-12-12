@testset "MorphSE test" begin
    println("------------------------------------------------")
    println("---------------- MorphSE Tests -----------------")
    # Dilate -- Start with a pixel in the middle and dilate in one go to fill up the full image
    n = rand(11:2:21) # choose random odd number
    mid = (n - 1) ÷ 2 + 1 # get median
    a = zeros(Int, n, n)
    a[mid, mid] = 1 # make 1 the pixel in the center
    se = IceFloeTracker.MorphSE.strel_box((n, n))
    @test IceFloeTracker.MorphSE.dilate(a, se) == ones(Int, n, n)

    # Bothat, opening, erode using output from Matlab
    A = zeros(Bool, 41, 41)
    A[(21 - 10):(21 + 10), (21 - 10):(21 + 10)] .= 1
    se = centered(IceFloeTracker.se_disk4())
    path = joinpath(test_data_dir, "morphSE")
    erode_withse_exp = readdlm(joinpath(path, "erode_withse1_exp.csv"), ',', Bool)
    bothat_withse_exp = readdlm(joinpath(path, "bothat_withse1_exp.csv"), ',', Bool)
    open_withse_exp = readdlm(joinpath(path, "open_withse1_exp.csv"), ',', Bool)
    @test open_withse_exp == IceFloeTracker.MorphSE.opening(A, se)
    @test erode_withse_exp == IceFloeTracker.MorphSE.erode(A, se)
    @test bothat_withse_exp == IceFloeTracker.MorphSE.bothat(A, se)
end
