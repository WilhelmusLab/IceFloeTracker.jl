@testset "Create Landmask" begin
    println("------------------------------------------------")
    println("------------ Create Landmask Test --------------")

    # define constants, maybe move to test config file
        landmask_file = "$(test_data_dir)/landmask.tiff"
        current_landmask_file = "$(test_data_dir)/current_landmask.png"
        lm_test_region = (1:800, 1:1500)
        matlab_landmask_file = "$(test_data_dir)/matlab_landmask.png"
        strel_file = "$(test_data_dir)/se.csv"
        struct_elem = BitArray(readdlm(strel_file, ',', Bool))
    
    # load images
        matlab_landmask = load(matlab_landmask_file)[lm_test_region...]
        lm_image = load(landmask_file)[lm_test_region...]
        test_image = load(truecolor_test_image_file)[lm_test_region...]

    # make landmask
        @time landmask = IceFloeTracker.create_landmask(lm_image, struct_elem)
    
    # test for percent difference in landmask images until ImageMorphology v0.5 is released
        @test 0.005 > sum(landmask .!== matlab_landmask) / prod(size(matlab_landmask))
        #@test (nothing === @test_approx_eq_sigma_eps landmask matlab_landmask [0, 0] .0001)

    # test for apply_landmask
        bw = ones(Int,3,3) # 3x3 box of ones
        mask = BitArray(rand([true, false],3,3)) # random 3x3 mask
        IceFloeTracker.apply_landmask!(bw,mask); # apply mask to bw in place
        @test .!mask == bw # check masked bw is the complement of mask
end
