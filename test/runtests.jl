using IceFloeTracker
using Images
using Test
using DelimitedFiles
using Dates

# Setting things up

## locate some files for the tests
test_data_dir = "./data"
test_image_file = "$(test_data_dir)/NE_Greenland_truecolor.2020162.aqua.250m.tiff"

## Get all test files filenames "test-*" in test folder and their corresponding names/label
alltests = [f for f in readdir() if startswith(f,"test-")]
testnames = [n[6:end-3] for n in alltests]

## Put the filenames to test below
to_test = 
        alltests # uncomment this line to run all tests or add individual files below 
        # [
        #  "test-create-landmask.jl",
        #  "test-create-cloudmask.jl",
        #  "test-normalize-image.jl",
        #  "test-persist.jl",
        #  "test-utils-padding.jl",
        #   ]

# Run the tests
@testset "IceFloeTracker.jl" begin
    for test in to_test
        include(test)
    end

    @testset "utils.jl pad utilities" begin
        println("-------------------------------------------------")
        println("--------------- Pad Image Tests -----------------")
        rad = 3
        val = 1
        simpleimg = collect(reshape(1:4,2,2))
        h,w = size(simpleimg)
        # 2×2 Matrix{Int64}:
        #  1  3
        #  2  4
        paddedimg = IceFloeTracker.add_padding(simpleimg, rad) # replicate with rad=3
        # 8×8 Matrix{Int64}:
        #  1  1  1  1  3  3  3  3
        #  1  1  1  1  3  3  3  3
        #  1  1  1  1  3  3  3  3
        #  1  1  1  1  3  3  3  3
        #  2  2  2  2  4  4  4  4
        #  2  2  2  2  4  4  4  4
        #  2  2  2  2  4  4  4  4
        #  2  2  2  2  4  4  4  4
        
        # test the corners
        @test paddedimg[1,1]     == 1 &&
              paddedimg[1,end]   == 3 &&
              paddedimg[end,1]   == 2 &&
              paddedimg[end,end] == 4 
        
        paddedimg = IceFloeTracker.add_padding(simpleimg, rad, :replicate) # again, replicate with rad = 3
        # test the corners
        @test paddedimg[1,1]     == 1 &&
              paddedimg[1,end]   == 3 &&
              paddedimg[end,1]   == 2 &&
              paddedimg[end,end] == 4 

        paddedimg = IceFloeTracker.add_padding(simpleimg, rad, :fill) # pad with zeros (by default) with rad=3
        # test the corners
        @test paddedimg[1,1]     == 
              paddedimg[1,end]   == 
              paddedimg[end,1]   == 
              paddedimg[end,end] == 0

        paddedimg = IceFloeTracker.add_padding(simpleimg, 3, :fill, 1) # pad with ones with rad=3
        @test paddedimg[1,1]     == 
              paddedimg[1,end]   == 
              paddedimg[end,1]   == 
              paddedimg[end,end] == 1
        
        # test sizing
        @test size(paddedimg) == (h+2*rad,w+2*rad)
        
        # test padding removal
        @test IceFloeTracker.remove_padding(paddedimg, rad) == simpleimg
    end    
end
