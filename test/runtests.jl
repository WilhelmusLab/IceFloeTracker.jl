using ArgParse: @add_arg_table!, ArgParseSettings, add_arg_group!, parse_args
using DataFrames
using Dates
using DelimitedFiles
using IceFloeTracker
using Images
using Random
using Test
using TestImages
using TiledIteration
using ZipFile
include("test_error_rate.jl")
include("config.jl")

function pad_string(str::String, total_length::Int=49, padding_char::Char='-')
    # Calculate the padding needed on each side
    left_padding = div(total_length - length(str), 2)
    right_padding = total_length - length(str) - left_padding

    # Pad the string
    padded_str = lpad(
        rpad(str, length(str) + right_padding, padding_char), total_length, padding_char
    )

    return padded_str
end

divline = "-"^49

macro ntestset(file, testblock)
    quote
        testset_name = basename($file)
        @testset "$testset_name" begin
            println(divline)
            println(pad_string(testset_name, length(divline)))
            $(esc(testblock))
        end
    end
end

# Setting things up (see config.jl)

## Get all test files filenames "test-*" in test folder and their corresponding names/label
alltests = [f for f in readdir() if startswith(f, "test-")]

## Put the filenames to test below

to_test = alltests # uncomment this line to run all tests or add individual files below
# to_test = [
# "test-branch.jl",
# "test-bridge.jl",
# "test-brighten.jl",
# "test-bwareamaxfilt.jl",
# "test-bwperim.jl",
# "test-bwtraceboundary.jl",
# "test-conditional-adaptive-histeq.jl",
# "test-create-cloudmask.jl",
# "test-create-landmask.jl",
# "test-crosscorr.jl",
# "test-discrim-ice-water.jl",
# "test-find-ice-labels.jl",
# "test-get-ice-masks.jl",
# "test-hbreak.jl",
# "test-imadjust.jl",
# "test-imcomplement.jl",
# "test-img-processing.jl",
# "test-latlon.jl",
# "test-long-tracker.jl",
# "test-matchcorr.jl",
# "test-misc.jl",
# "test-morph-fill.jl",
# "test-morphSE.jl",
# "test-normalize-image.jl",
# "test-persist.jl",
# "test-psi-s.jl",
# "test-reconstruction-workflows.jl",
# "test-regionprops.jl",
# "test-regionprops-labeled.jl",
# "test-register.jl",
# "test-register-utils.jl",
# "test-regularize-final.jl",
# "test-resample-boundary.jl",
# "test-rotation.jl",
# "test-segmentation-a.jl",
# "test-segmentation-b.jl",
# "test-segmentation-f.jl",
# "test-segmentation-lopez-acosta-2019.jl",
# "test-segmentation-lopez-acosta-2019-tiling.jl",
# "test-segmentation-watershed.jl",
# "test-special-strels.jl",
# "test-tiled-ice-labels.jl",
# "test-tilingutils.jl",
# "test-utils-imextendedmin.jl",
# "test-utils-padding.jl",
# "test-watershed-workflows.jl",
# ]

# Run the tests
@testset verbose = true "IceFloeTracker.jl" begin
    for test in to_test
        include(test)
    end
end
