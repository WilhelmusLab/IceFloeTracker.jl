"""
Module for command-line interface for ice floe tracking.
"""
module CLI
using ..Segmentation
using Images

function (p::IceFloeSegmentationAlgorithm)(
    outpath::AbstractString,
    args::AbstractString...;
    intermediates_directory::Union{Nothing,AbstractString}=nothing,
    intermediates_targets::Vector{<:AbstractString}=AbstractString[],
)
    intermediate_results_callback = intermediate_results_saver_factory(
        intermediates_directory, intermediates_targets
    )
    loaded_inputs = []
    for name in args
        input = load(name)
        if occursin("landmask", name)
            input = load(name) .|> Gray .|> (x -> x .> 0.1) .|> Gray
        else
            input = load(name)
        end
        @info "Loaded input $(name): $(size(input)) $(input[1])..."
        push!(loaded_inputs, input)
    end
    output = p(loaded_inputs...; intermediate_results_callback)

    save(outpath, reinterpret(Gray{N0f64}, labels_map(output)))
    return output
end

function intermediate_results_saver_factory(
    directory::Union{Nothing,AbstractString}, targets::Vector{<:AbstractString}
)
    if isnothing(directory)
        return (; kwargs...) -> nothing
    end

    function save_intermediate_results(; kwargs...)
        for name in targets
            root = Symbol(splitext(name)[1])
            if haskey(kwargs, root)
                mkpath(directory)
                save(joinpath(directory, name), kwargs[root])
            end
        end
    end

    return save_intermediate_results
end

end # module CLI