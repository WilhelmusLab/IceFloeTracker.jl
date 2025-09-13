module NotebookToDocumenter
# Based on https://github.com/marius311/CMBLensing.jl/blob/v0.10.1/docs/make.jl
# TODO:
# - replace convert_equations! with jinja template function

function notebooks_to_documenter_md(directory)
    for path in readdir_recursive(directory)
        if endswith(path, ".ipynb")
            @info "Converting $path to markdown"
            notebook_to_documenter_md(path)
        end
    end
end

function notebook_to_documenter_md(file)
    new_file = file |> convert_to_markdown |> convert_equations!
    return new_file
end

function convert_to_markdown(file)
    template_path = joinpath(dirname(@__FILE__), "documenter.tpl")
    run(
        `pipx run --spec nbconvert jupyter-nbconvert $file --to markdown --template $template_path`,
    )
    new_file = replace(file, "ipynb" => "md")
    return new_file
end

function convert_equations!(file)
    contents = read(file, String)
    contents = replace(contents, r"\$\$(.*?)\$\$"s => s"""```math
    \g<1>
    ```""")
    contents = replace(contents, r"\* \$(.*?)\$" => s"* ``\g<1>``") # starting a line with inline math screws up tex2jax for some reason
    write(file, contents)
    return file
end

function readdir_recursive(directory)
    result::Vector{String} = String[]
    for (root, _, files) in walkdir(directory)
        for file in files
            push!(result, joinpath(root, file))
        end
    end
    return result
end

end

run(`rsync --recursive --delete docs/src/ docs/prebuild/`)
NotebookToDocumenter.notebooks_to_documenter_md("docs/prebuild")

using Documenter
using IceFloeTracker

makedocs(;
    sitename="IceFloeTracker.jl",
    format=Documenter.HTML(; size_threshold=nothing),
    modules=[IceFloeTracker],
    doctest=false,
    warnonly=true,
    source="prebuild/",
)

deploydocs(;
    repo="github.com/WilhelmusLab/IceFloeTracker.jl.git",
    push_preview=true,
    versions=nothing,
)
