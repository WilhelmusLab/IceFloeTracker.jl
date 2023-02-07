#!/usr/bin/env julia
using Pkg
Pkg.activate(@__DIR__)

using ArgParse
using IceFloeTracker

function main(args)
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "fetchdata"
        help = "Fetch source data for ice floe tracking"
        action = :command

        "landmask"
        help = "Generate land mask images"
        action = :command

        "cloudmask"
        help = "Generate cloud mask images"
        action = :command

        "--extract_features", "-e"
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

    @add_arg_table! settings["extract_features"] begin
        "--input", "-i"
        help = "Input image directory"
        required = true

        "--output", "-o"
        help = "Output image directory"
        required = true

        "--area_threshold", "-a"
        help = "Minimum and maximum area of ice floes to extract"
        required = true

        "--features", "-f"
        help = "Features to extract"
        required = false
        default = ["centroid", "area", "major_axis_length", "minor_axis_length", "convex_area", "bbox"]
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
    add_arg_table!(settings["cloudmask"], landmask_cloudmask_args...)

    parsed_args = parse_args(args, settings; as_symbols=true)

    command = parsed_args[:_COMMAND_]
    command_args = parsed_args[command]
    command_func = getfield(IceFloeTracker, Symbol(command))
    command_func(; command_args...)
    return nothing
end

main(ARGS)
