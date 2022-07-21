using IceFloeTracker
using Images
using Test
using DelimitedFiles
using Dates

# Setting things up

## locate some files for the tests
test_data_dir = "./test_inputs"
test_output_dir = "./test_outputs"
truecolor_test_image_file = "$(test_data_dir)/NE_Greenland_truecolor.2020162.aqua.250m.tiff"
reflectance_test_image_file = "$(test_data_dir)/NE_Greenland.2020162.aqua.250m.tiff"
landmask_file = "$(test_data_dir)/landmask.tiff"
current_landmask_file = "$(test_data_dir)/current_landmask.png"
test_region = (1:2707, 1:4458)

## Get all test files filenames "test-*" in test folder and their corresponding names/label
alltests = [f for f in readdir() if startswith(f, "test-")]
testnames = [n[6:(end - 3)] for n in alltests]

## Put the filenames to test below
to_test = #alltests # uncomment this line to run all tests or add individual files below 
[
    "test-create-landmask.jl",
    "test-create-cloudmask.jl",
    "test-normalize-image.jl",
    #    "test-persist.jl",
    #    "test-utils-padding.jl",
    #"test-discrim-ice-water.jl",
]

# Run the tests
@testset "IceFloeTracker.jl" begin
    for test in to_test
        include(test)
    end
end
