@testset "MorphSE test" begin
    println("------------------------------------------------")
    println("---------------- MorphSE Tests -----------------")
    # Dilate -- Start with a pixel in the middle and dilate in one go to fill up the full image
    n = rand(11:2:21) # choose random odd number
    mid = (n - 1) ÷ 2 + 1 # get median
    a = zeros(Int, n, n)
    a[mid, mid] = 1 # make 1 the pixel in the center
    se = IceFloeTracker.strel_box((n, n))
    @test IceFloeTracker.dilate(a, se) == ones(Int, n, n)

    # Bothat, opening, erode, filling holes, reconstruction using output from Matlab
    A = zeros(Bool, 41, 41)
    A[(21 - 10):(21 + 10), (21 - 10):(21 + 10)] .= 1
    I = falses(8, 8)
    I[1:8, 3:6] .= 1
    I[[CartesianIndex(4, 4), CartesianIndex(5, 5)]] .= 0
    I
    se = centered(IceFloeTracker.se_disk4())

    #read in expected files from MATLAB
    path = joinpath(test_data_dir, "morphSE")
    erode_withse_exp = readdlm(joinpath(path, "erode_withse1_exp.csv"), ',', Bool)
    bothat_withse_exp = readdlm(joinpath(path, "bothat_withse1_exp.csv"), ',', Bool)
    open_withse_exp = readdlm(joinpath(path, "open_withse1_exp.csv"), ',', Bool)
    reconstruct_exp = readdlm(joinpath(path, "reconstruct_exp.csv"), ',', Int64)
    matrix_A = readdlm(joinpath(path, "mat_a.csv"), ',', Int64)
    matrix_B = readdlm(joinpath(path, "mat_b.csv"), ',', Int64)
    filled_holes_exp = readdlm(joinpath(path, "filled_holes.csv"), ',', Int64)

    #run tests
    @test open_withse_exp == IceFloeTracker.opening(A, se)
    @test erode_withse_exp == IceFloeTracker.erode(A, se)
    @test bothat_withse_exp == IceFloeTracker.bothat(A, se)
    @test reconstruct_exp ==
        IceFloeTracker.mreconstruct(IceFloeTracker.dilate, matrix_B, matrix_A)
    @test filled_holes_exp == IceFloeTracker.fill_holes(I)
end
