"""
    run_segmentation_over_multiple_cases(
        data_loader::ValidationDataLoader,
        case_filter::Function,
        algorithm::IceFloeSegmentationAlgorithm;
        output_directory::Union{AbstractString,Nothing}=nothing,
        result_images_to_save::Union{AbstractArray{Symbol},Nothing}=nothing,
    )

Run the `algorithm::IceFloeSegmentationAlgorithm` over each of the `cases` and return a DataFrame of the results.

Inputs:
- `data_loader`: a ValidationDataLoader to load the data from a validated data set
- `case_filter`: a function which returns a boolean when given the data_loader's metadata values
  and determines which cases are included (true) or excluded (false) from processing
- `algorithm`: an instantiated `IceFloeSegmentationAlgorithm` which will be called on each case
- `output_directory`: optional – path to save intermediate and final outputs
- `result_images_to_save`: optional – symbols of the intermediate results from `algorithm` which should be saved

Returns:
- A DataFrame with the results including a :success boolean, any :error messages, and the original metadata.

"""

function run_segmentation_over_multiple_cases(
    data_loader::ValidationDataLoader,
    case_filter::Function,
    algorithm::IceFloeSegmentationAlgorithm;
    output_directory::Union{AbstractString,Nothing}=nothing
)::DataFrame
    dataset = data_loader(; case_filter)
    @info dataset.metadata
    results = []
    for case::ValidationDataCase in dataset
        let name, datestamp, validated, measured, success, error, comparison
            name = case.name
            datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
            validated = case.validated_labeled_floes

            if !isnothing(output_directory)
                intermediate_results_callback = save_results_callback(
                    joinpath(
                        output_directory,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)",
                    );
                )
            else
                intermediate_results_callback = nothing
            end

            @info "starting $(name)"
            try
                measured = algorithm(
                    RGB.(case.modis_truecolor),
                    RGB.(case.modis_falsecolor),
                    case.modis_landmask;
                    intermediate_results_callback,
                )
                @info "$(name) succeeded"
                success = true
                error = nothing
            catch error
                @warn "$(name) failed: $(error)"
                success = false
                measured = nothing
            end

            comparison = segmentation_comparison(; validated, measured)

            # Store the aggregate results
            push!(
                results,
                merge((; name, success, error), comparison, NamedTuple(case.metadata)),
            )
            if !isnothing(intermediate_results_callback) && !isnothing(validated)
                intermediate_results_callback(;
                    segment_mean_truecolor_validated=map(
                        i -> segment_mean(validated, i), labels_map(validated)
                    ),
                )
            end
        end
    end
    results_df = DataFrame(results)
    @info results_df
    return results_df
end

"""
    save_results_callback(
        path;
        extension,
        names::Union{AbstractArray{Symbol},Nothing}
    )::Function

Returns a function which saves any images passed into it as keyword arguments.

# Example
```julia-repl
julia> callback = save_results_callback("/tmp/path/to/directory")
julia> image = Gray.([1 1 0 0 1 0 1])
julia> callback(;image_name=image)
```
... saves `image` to `/tmp/path/to/directory/image_name.png`.
"""
function save_results_callback(
    directory::AbstractString;
    extension::AbstractString=".png",
    names::Union{AbstractArray{Symbol},Nothing}=nothing,
)
    function callback(; kwargs...)
        mkpath(directory)
        for (name, image) in kwargs
            (names === nothing || name ∈ names) || continue
            path = joinpath(directory, String(name) * extension)
            if typeof(image) <: AbstractArray{Bool}
                image = Gray.(image)
            end
            if typeof(image) <: AbstractArray{<:Colorant}
                try
                    save(path, image)
                catch e
                    @warn "an unexpected error occured saving $name: $e"
                end
            else
                @debug "skipping $(name) – not an image we can save"
            end
        end
    end
    return callback
end