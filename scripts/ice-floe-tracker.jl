#!/usr/bin/env julia
using Pkg
# Pkg.activate(@__DIR__)

using ArgParse
using IceFloeTracker

function main(args)
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "fetchdata"
        help = "Fetch source data for ice floe tracking"
        action = :command

        "landmask"
        help = "generates land mask images"
        action = :command

        "prep"
        help = "Preprocess input images"
        action = :command

        "seg"
        help = "Do segmentation on preprocessed images"
        action = :command

        "fext"
        help = "Do feature extraction on segmented images"
        action = :command

        "track"
        help = "Do tracking procedure using extracted features from segmented images"
        action = :command

        "ltrack"
        help = "Do long tracking procedure"
        action = :command
    end

    @add_arg_table! settings["fetchdata"] begin
        "output"
        help = "output image directory"
        required = true
    end

    @add_arg_table! settings["landmask"] begin
        "output"
        help = "output image directory"
        required = true
    end

    # metadata requirements might change later
    command_common_args = [
        "metadata",
        Dict(:help => "image metadata file", :required => false),
        "input",
        Dict(:help => "input image directory", :required => true),
        "output",
        Dict(:help => "output image directory", :required => true),
    ]
    # add_arg_table!(settings["landmask"], command_common_args...)
    add_arg_table!(settings["prep"], command_common_args...)
    add_arg_table!(settings["seg"], command_common_args...)
    add_arg_table!(settings["fext"], command_common_args...)
    add_arg_table!(settings["track"], command_common_args...)

    parsed_args = parse_args(args, settings; as_symbols=true)

    command = parsed_args[:_COMMAND_]
    command_args = parsed_args[command]
    command_func = getfield(IceFloeTracker, Symbol(parsed_args[:_COMMAND_]))

    command_func(; command_args...)
    return nothing
end

main(ARGS)
