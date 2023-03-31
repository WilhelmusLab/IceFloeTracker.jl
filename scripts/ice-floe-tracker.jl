#!/usr/bin/env julia
using Pkg
Pkg.activate(@__DIR__)

using ArgParse
using IceFloeTracker

function main(args)
    settings = ArgParseSettings(; autofix_names=true)

    @add_arg_table! settings begin
        "fetchdata"
        help = "Fetch source data for ice floe tracking"
        action = :command

        "landmask"
        help = "Generate land mask images"
        action = :command

        "preprocess"
        help = "Preprocess truecolor/reflectance images"
        action = :command

        "extractfeatures"
        help = "Extract ice floe features from segmented floe image"
        action = :command
    end

    @add_arg_table! settings["fetchdata"] begin
        "output"
        help = "Output image directory"
        required = true
    end

    landmask_cloudmask_args = [
        "input",
        Dict(:help => "Input image directory", :required => true),
        "output",
        Dict(:help => "Output image directory", :required => true),
    ]

    @add_arg_table! settings["preprocess"] begin
        "--truedir", "-t"
        help = "Truecolor image directory"
        required = true

        "--refdir", "-r"
        help = "Reflectance image directory"
        required = true

        "--lmdir", "-l"
        help = "Land mask image directory"
        required = true

        "--output", "-o"
        help = "Output directory"
        required = true
        end

    @add_arg_table! settings["extractfeatures"] begin
        "--input", "-i"
        help = "Input image directory"
        required = true

        "--output", "-o"
        help = "Output image directory"
        required = true

        "--minarea"
        help = "Minimum area (in pixels) of ice floes to extract"
        required = false
        default = "300"

        "--maxarea"
        help = "Maximum area (in pixels) of ice floes to extract"
        required = false
        default = "90000"

        "--features", "-f"
        help = """Features to extract. Format: "feature1 feature2". For an extensive list of extractable features see https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.regionprops:~:text=The%20following%20properties%20can%20be%20accessed%20as%20attributes%20or%20keys"""
        required = false
        default = "centroid area major_axis_length minor_axis_length convex_area bbox orientation"
    end

    command_common_args = [
        "metadata",
        Dict(:help => "Image metadata file", :required => true),
        "input",
        Dict(:help => "Input image directory", :required => true),
        "output",
        Dict(:help => "Output image directory", :required => true),
    ]

    add_arg_table!(settings["landmask"], landmask_cloudmask_args...)
    
    parsed_args = parse_args(args, settings; as_symbols=true)

    command = parsed_args[:_COMMAND_]
    command_args = parsed_args[command]
    command_func = getfield(IceFloeTracker.Pipeline, Symbol(command))
    command_func(; command_args...)
    return nothing
end

main(ARGS)
