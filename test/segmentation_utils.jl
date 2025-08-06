"""
Run the `algorithm::IceFloeSegmentationAlgorithm` over each of the `cases` and return a DataFrame of the results.
- Each `case` should be a ValidationDataCase.
- If `output_directory` is defined, then save the output segmentations and images to the directory. 
"""
function run_segmentation_over_multiple_cases(
    cases,
    algorithm::IceFloeSegmentationAlgorithm;
    output_directory::Union{AbstractString,Nothing}=nothing,
)::DataFrame
    results = []
    for case::ValidationDataCase in cases
        validated = case.validated_labeled_floes
        name = case.name
        let measured, success, error
            @info "starting $(name)"
            try
                measured = algorithm(
                    RGB.(case.modis_truecolor),
                    RGB.(case.modis_falsecolor),
                    case.modis_landmask,
                )
                @info "$(name) succeeded"
                success = true
                error = nothing
            catch error
                @warn "$(name) failed: $(error)"
                success = false
                measured = nothing
            end

            # Store the aggregate results
            push!(
                results,
                merge(
                    (; name, success, error), segmentation_comparison(; validated, measured)
                ),
            )

            if !isnothing(output_directory)
                mkpath(output_directory)
                datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
                !isnothing(measured) && save(
                    joinpath(
                        output_directory,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)-mean-labels.png",
                    ),
                    map(i -> segment_mean(measured, i), labels_map(measured)),
                )
                !isnothing(validated) && save(
                    joinpath(
                        output_directory,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)-validated-mean-labels.png",
                    ),
                    map(i -> segment_mean(validated, i), labels_map(validated)),
                )
                !isnothing(case.modis_truecolor) && save(
                    joinpath(
                        output_directory,
                        "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)-truecolor.png",
                    ),
                    case.modis_truecolor,
                )
            end
        end
    end
    results_df = DataFrame(results)
    return results_df
end
