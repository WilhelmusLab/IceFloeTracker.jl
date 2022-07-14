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
tests = Dict(testfile=>name for (testfile,name) in zip(alltests,testnames))

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
end

