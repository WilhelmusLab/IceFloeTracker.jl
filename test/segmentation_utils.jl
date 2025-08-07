"""
Run the `algorithm::IceFloeSegmentationAlgorithm` over each of the `cases` and return a DataFrame of the results.
- Each `case` should be a ValidationDataCase.
- If `output_directory` is defined, then save the output segmentations and images to the directory. 
"""
function run_segmentation_over_multiple_cases(
    cases,
    algorithm::IceFloeSegmentationAlgorithm;
    output_directory::Union{AbstractString,Nothing}=nothing,
    result_images_to_save::Union{AbstractArray{Symbol},Nothing}=nothing,
)::DataFrame
    results = []
    for case::ValidationDataCase in cases
        datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")

        let name, validated, measured, success, error, comparison
            name = case.name
            validated = case.validated_labeled_floes

            if !isnothing(output_directory)
                intermediate_results_callback = save_results_callback(
                    joinpath(
                        output_directory,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)",
                    );
                    names=result_images_to_save,
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
            push!(results, merge((; name, success, error), comparison))

            # Save the validated results we have them
            if !isnothing(validated) && !isnothing(intermediate_results_callback)
                intermediate_results_callback(;
                    validated_mean_labels=map(
                        i -> segment_mean(validated, i), labels_map(validated)
                    ),
                )
            end
        end
    end
    results_df = DataFrame(results)
    return results_df
end
