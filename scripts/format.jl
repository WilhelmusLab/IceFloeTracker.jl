#!/usr/bin/env julia
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Git
using JuliaFormatter

function main()
    staged = readlines(`$git diff --cached --name-only --diff-filter=AM`)
    staged = [f for f in staged if endswith(f, "jl")]
    format(staged, BlueStyle())
    status = (format(staged, BlueStyle();)) ? 0 : 1
    return exit(status)
end

main()
