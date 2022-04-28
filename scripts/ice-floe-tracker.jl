#!/usr/bin/env julia --project=@.

using ArgParse
using IceFloeTracker

function main(args)
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "fetchdata"
        help = "Fetches source data for ice floe tracking"
        action = :command

        "landmask"
        help = "generates land mask images"
        action = :command

        "cloudmask"
        help = "generates cloud mask images"
        action = :command
    end

    @add_arg_table! settings["fetchdata"] begin
        "output"
        help = "output image directory"
        required = true
    end

    command_common_args = [
        "metadata",
        Dict(:help => "image metadata file", :required => true),
        "input",
        Dict(:help => "input image directory", :required => true),
        "output",
        Dict(:help => "output image directory", :required => true),
    ]
    add_arg_table!(settings["landmask"], command_common_args...)
    add_arg_table!(settings["cloudmask"], command_common_args...)

    parsed_args = parse_args(args, settings; as_symbols=true)

    command = parsed_args[:_COMMAND_]
    command_args = parsed_args[command]
    command_func = getfield(IceFloeTracker, Symbol(parsed_args[:_COMMAND_]))

    command_func(; command_args...)
    return nothing
end

main(ARGS)
