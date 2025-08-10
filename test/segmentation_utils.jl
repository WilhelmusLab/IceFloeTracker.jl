"""
Compares two SegmentedImages and returns values describing how similar the segmentations are.

Treats the segment labeled `0` as background.

Measures:
- normalized_{validated,measured}_area: fraction of image covered by segments
- fractional_intersection: fraction of the validated segments covered by measured segments
"""
function segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}=nothing,
    measured::Union{SegmentedImage,Nothing}=nothing,
)::NamedTuple
    if !isnothing(validated)
        validated_binary = Gray.(labels_map(validated) .> 0)
        validated_area = sum(channelview(validated_binary))
        normalized_validated_area = validated_area / length(channelview(validated_binary))
    else
        normalized_validated_area = missing
    end

    if !isnothing(measured)
        measured_binary = Gray.(labels_map(measured) .> 0)
        measured_area = sum(channelview(measured_binary))
        normalized_measured_area = measured_area / length(channelview(measured_binary))
    else
        normalized_measured_area = missing
    end

    if !isnothing(validated) && !isnothing(measured)
        intersection = Gray.(channelview(measured_binary) .&& channelview(validated_binary))
        fractional_intersection =
            fractional_intersection = sum(channelview(intersection)) / validated_area
    else
        fractional_intersection = missing
    end

    return (; normalized_validated_area, normalized_measured_area, fractional_intersection)
end

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

function run_segmentation_over_multiple_cases(
    data_loader::ValidationDataLoader,
    case_filter::Function,
    algorithm::IceFloeSegmentationAlgorithm;
    output_directory::Union{AbstractString,Nothing}=nothing,
)::@NamedTuple{metadata::DataFrame, results::DataFrame}
    dataset = data_loader(; case_filter)
    @info dataset.metadata
    results = run_segmentation_over_multiple_cases(
        dataset.data, algorithm; output_directory
    )
    @info results
    return (; metadata=dataset.metadata, results)
end

function test_all_cases_ran_without_crashing(
    outputs::@NamedTuple{metadata::DataFrame, results::DataFrame};
    success_column::Symbol=:success,
)
    # All cases from the metadata are included in the results
    @test nrow(outputs.results) == nrow(outputs.metadata)

    # ... and each of them succeeded
    successes = subset(outputs.results, success_column => ByRow(==(true)))
    @test nrow(outputs.results) == nrow(successes)
end