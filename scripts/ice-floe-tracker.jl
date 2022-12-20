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
    end

    @add_arg_table! settings["fetchdata"] begin
        "output"
        help = "Output image directory"
        required = true
    end

    landmask_args = [
        "input",
        Dict(:help => "Input image directory", :required => true),
        "output",
        Dict(:help => "Output image directory", :required => true),
    ]

    command_common_args = [
        "metadata",
        Dict(:help => "Image metadata file", :required => true),
        "input",
        Dict(:help => "Input image directory", :required => true),
        "output",
        Dict(:help => "Output image directory", :required => true),
    ]

    add_arg_table!(settings["landmask"], landmask_args...)
    add_arg_table!(settings["cloudmask"], command_common_args...)

    parsed_args = parse_args(args, settings; as_symbols=true)

    command = parsed_args[:_COMMAND_]
    command_args = parsed_args[command]
    command_func = getfield(IceFloeTracker, Symbol(command))
    command_func(; command_args...)
    return nothing
end

main(ARGS)
