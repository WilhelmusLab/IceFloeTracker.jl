using Documenter
using IceFloeTracker

# Based on https://github.com/marius311/CMBLensing.jl/blob/v0.10.1/docs/make.jl

function convert_notebooks(directory; converter)
    for path in readdir_recursive(directory)
        if endswith(path, ".ipynb")
            @info "Converting $path"
            converter(path)
        end
    end
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

function convert_example_blocks!(file)
    contents = read(file, String)
    contents = replace(contents, r"```julia(.*?)```"s => s"""```@example _page-environment
    \g<1>
    ```""")
    write(file, contents)
    return file
end

function add_colab_link!(; kwargs...)
    return (file -> add_colab_link!(file; kwargs...))
end

function add_colab_link!(
    file;
    path_resolver="docs/prebuild/" => "docs/src/",
    username="username",
    repo="repo",
    ref="main",
    extension=".ipynb",
    colab_badge_url="https://colab.research.google.com/assets/colab-badge.svg",
    alt_text="Open this notebook in Colab",
)
    source_relative_path = replace(file, path_resolver)
    target_relative_path = splitext(source_relative_path)[1] * extension
    colab_link = joinpath(
        "https://colab.research.google.com/github/",
        username,
        repo,
        "blob",
        ref,
        target_relative_path,
    )
    colab_badge = "[![$alt_text]($colab_badge_url)]($colab_link)"
    contents = read(file, String)
    contents = colab_badge * "\n" * contents
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

username = "WilhelmusLab"
repo = "IceFloeTracker.jl"
ref = get(ENV, "GITHUB_REF", "main")

run(`rsync --recursive --delete docs/src/ docs/prebuild/`)
convert_notebooks(
    "docs/prebuild";
    converter=file ->
        file |>
        convert_to_markdown |>
        convert_equations! |>
        convert_example_blocks! |>
        add_colab_link!(; username, repo, ref),
)

makedocs(;
    sitename="IceFloeTracker.jl",
    format=Documenter.HTML(; size_threshold=nothing),
    modules=[IceFloeTracker],
    doctest=false,
    warnonly=true,
    source="prebuild/",
    pages=[
        "IceFloeTracker.jl" => "index.md",
        "preprocessing.md",
        "segmentation.md",
        "tracking.md",
        "Tutorials" => [
            "tutorials/lopez-acosta-2019-workflow.md",
            "tutorials/preprocessing-workflow.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/WilhelmusLab/IceFloeTracker.jl.git",
    devbranch="main",
    push_preview=true,
    versions=nothing,
)
