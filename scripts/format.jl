#!/usr/bin/env julia --project=@.
using ArgParse
using JuliaFormatter

function check_format(filepath::AbstractString)::Bool
    open(filepath, "r") do file
        a = read(file, String)
        b = format_text(a)

        return (a == b) ? true : false
    end
end

function format(filepath::AbstractString)::Bool
    return format_file(filepath)
end

function main(args)
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "--check", "-c"
        help = "Check to see if formatting is OK"
        action = "store_true"
    end

    parsed_args = parse_args(args, settings; as_symbols = true)

    julia_files = []
    for (root, dirs, files) in walkdir(pwd())
        for file in files
            if endswith(file, ".jl")
                push!(julia_files, joinpath(root, file))
            end
        end
    end

    status = 0
    for julia_file in julia_files
        if parsed_args[:check]
            if !check_format(julia_file)
                println("$(julia_file) is poorly formatted")
                status = 1
            end
        else
            if !format(julia_file)
                println("$(julia_file) formatted")
                status = 1
            end
        end
    end

    exit(status)
end

main(ARGS)
