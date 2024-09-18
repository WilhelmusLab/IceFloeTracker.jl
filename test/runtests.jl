using ArgParse: @add_arg_table!, ArgParseSettings, add_arg_group!, parse_args
using DataFrames
using Dates
using DelimitedFiles
using IceFloeTracker
using ImageTransformations: imrotate
using Images
using Random
using Test
using TestImages
using TiledIteration
include("test_error_rate.jl")
include("config.jl")

# Setting things up (see config.jl)

## Get all test files filenames "test-*" in test folder and their corresponding names/label
alltests = [f for f in readdir() if startswith(f, "test-")]
testnames = [n[6:(end - 3)] for n in alltests]

## Put the filenames to test below

to_test = #alltests # uncomment this line to run all tests or add individual files below
[
    # "test-create-landmask.jl",
    # "test-create-cloudmask.jl",
    # "test-normalize-image.jl",
    # "test-persist.jl",
    # "test-utils-padding.jl",
    # "test-discrim-ice-water.jl",
    # "test-find-ice-labels.jl",
    # "test-segmentation-a.jl",
    # "test-segmentation-b.jl",
    # "test-segmentation-watershed.jl",
    # "test-segmentation-f.jl",
    # "test-bwtraceboundary.jl",
    # "test-resample-boundary.jl",
    # "test-regionprops.jl",
    # "test-psi-s.jl",
    # "test-crosscorr.jl"
    # "test-bwperim.jl",
    # "test-bwareamaxfilt.jl"
    # "test-register-mismatch.jl",
    # "test-utils-imextendedmin.jl",
    # "test-morphSE.jl",
    # "test-hbreak.jl",
    # "test-bridge.jl",
    # "test-branch.jl"
    # "test-pipeline.jl"
]

# Run the tests
@testset verbose = true "IceFloeTracker.jl" begin
    for test in to_test
        include(test)
    end
end
