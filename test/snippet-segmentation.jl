@testsnippet Segmentation begin
    using Images: labels_map, segment_mean, Colorant, save
    using DataFrames: DataFrame

    """
        run_and_validate_segmentation(
            dataset::ValidationDataSet,
            algorithm::IceFloeSegmentationAlgorithm;
            output_directory::Union{AbstractString,Nothing}=nothing,
        )

    Run the `algorithm::IceFloeSegmentationAlgorithm` over each of the cases in `dataset` and return a DataFrame of the results.

    Inputs:
    - `dataset`: a ValidationDataSet
    - `algorithm`: an instantiated `IceFloeSegmentationAlgorithm` which will be called on each case
    - `output_directory`: optional – path to save intermediate and final outputs

    Returns:
    - A DataFrame with the results including a :success boolean, any :error messages, and the original metadata.

    """
    function run_and_validate_segmentation(
        dataset::ValidationDataSet,
        algorithm::IceFloeSegmentationAlgorithm;
        output_directory::Union{AbstractString,Nothing}=nothing,
    )::DataFrame
        @info dataset.metadata
        results = []
        for case::ValidationDataCase in dataset
            @info "starting $(case.name)"
            results_row = run_and_validate_segmentation(case, algorithm; output_directory)
            push!(results, results_row)
        end
        results_df = DataFrame(results)
        @info results_df
        return results_df
    end

    """
        run_and_validate_segmentation(
            case::ValidationDataCase,
            algorithm::IceFloeSegmentationAlgorithm;
            output_directory::Union{AbstractString,Nothing}=nothing,
        )

    Run the `algorithm::IceFloeSegmentationAlgorithm` on the `case` and return a NamedTuple of the validation results.

    - `case`: ValidationDataCase to be processed by the algorithm
    - `algorithm`: an instantiated `IceFloeSegmentationAlgorithm` which will be called on each case
    - `output_directory`: optional – path to save intermediate and final outputs


    Results include:
    - `name` – name of the `case`
    - `success` – whether the algorithm ran without throwing an error
    - `error` – if there was an error, what the error message was
    - outputs from `segmentation_summary` including `labeled_fraction`
    - outputs from `segmentation_comparison` of the measured and validated dataset, including `precision`, `recall` and `F_score`.


    """
    function run_and_validate_segmentation(
        case::ValidationDataCase,
        algorithm::IceFloeSegmentationAlgorithm;
        output_directory::Union{AbstractString,Nothing}=nothing,
    )
        let name, datestamp, validated, measured, success, error, comparison
            name = case.name
            validated = case.validated_labeled_floes
            if !isnothing(output_directory)
                intermediate_results_callback = save_results_callback(
                    output_directory, case, algorithm;
                )
            else
                intermediate_results_callback = nothing
            end
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
                @error "An error occurred." exception = (error, catch_backtrace())
                success = false
                measured = nothing
            end
            summary = segmentation_summary(measured)
            comparison = segmentation_comparison(; validated, measured)
            results = merge(
                (; name, success, error), comparison, summary, NamedTuple(case.metadata)
            )
            if !isnothing(intermediate_results_callback) && !isnothing(validated)
                intermediate_results_callback(;
                    segment_mean_truecolor_validated=map(
                        i -> segment_mean(validated, i), labels_map(validated)
                    ),
                )
            end
            return results
        end
    end

    """
        results_invariant_for(
            target_type::Union{Function,Type},
            baseline::NamedTuple,
            algorithm::IceFloeSegmentationAlgorithm,
            case::ValidationDataCase,
        )::Bool
        
    Runs `algorithm` on `case` using `target_type` to cast images; returns true if results are within 1% of the `baseline`.


    """
    function results_invariant_for(
        target_type::Union{Function,Type}...;
        baseline::NamedTuple,
        algorithm::IceFloeSegmentationAlgorithm,
        case::ValidationDataCase,
        output_directory::AbstractString="./test_outputs",
        rtol::Real=0.01,
    )::Bool
        intermediate_results_callback = save_results_callback(
            joinpath(
                output_directory,
                "segmentation-$(typeof(algorithm))-$(join(target_type,"∘"))-$(Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS"))",
            );
        )
        casting_function = ∘(target_type...)
        segments = algorithm(
            casting_function.(case.modis_truecolor),
            casting_function.(case.modis_falsecolor),
            casting_function.(case.modis_landmask);
            intermediate_results_callback,
        )
        (; segment_count, labeled_fraction) = segmentation_summary(segments)
        (; recall, precision, F_score) = segmentation_comparison(
            case.validated_labeled_floes, segments
        )

        segment_count_pass = ≈(segment_count, baseline.segment_count; rtol)
        labeled_fraction_pass = ≈(labeled_fraction, baseline.labeled_fraction; rtol)
        recall_pass = ≈(recall, baseline.recall; rtol)
        precision_pass = ≈(precision, baseline.precision; rtol)
        F_score_pass = ≈(F_score, baseline.F_score; rtol)

        result = all([
            segment_count_pass,
            labeled_fraction_pass,
            recall_pass,
            precision_pass,
            F_score_pass,
        ])
        if !result
            @show segments
            @info "$(join(target_type,"∘")) failed"
            !segment_count_pass && @show (segment_count, baseline.segment_count)
            !labeled_fraction_pass && @show (labeled_fraction, baseline.labeled_fraction)
            !recall_pass && @show (recall, baseline.recall)
            !precision_pass && @show (precision, baseline.precision)
            !F_score_pass && @show (F_score, baseline.F_score)
        end

        return result
    end

    """
        save_results_callback(
            directory::AbstractString,
            case::ValidationDataCase,
            algorithm::IceFloeSegmentationAlgorithm;
            extension::AbstractString=".png",
        )::Function

    Returns a function which saves any images which are passed into it as keyword arguments.
    Creates a subdirectory based on the current time, the `case` and `algorithm`.

    Inputs:
    - `directory`: base directory where images will be stored
    - `case`: ValidationDataCase with metadata which are used to name a subdirectory
    - `algorithm`: IceFloeSegmentationAlgorithm which is used in the subdirectory name.
    """
    function save_results_callback(
        directory::AbstractString,
        case::ValidationDataCase,
        algorithm::IceFloeSegmentationAlgorithm;
        extension::AbstractString=".png",
    )::Function
        datestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
        name = case.name
        return save_results_callback(
            joinpath(directory, "segmentation-$(typeof(algorithm))-$(name)-$(datestamp)");
            extension,
        )
    end

    """
        save_results_callback(
            path;
            extension,
        )::Function

    Returns a function which saves any images which are passed into it as keyword arguments.

    # Example
    ```julia-repl
    julia> callback = save_results_callback("/tmp/path/to/directory")
    julia> image = Gray.([1 1 0 0 1 0 1])
    julia> callback(;image_name=image)
    ```
    ... saves `image` to `/tmp/path/to/directory/image_name.png`.
    """
    function save_results_callback(
        directory::AbstractString; extension::AbstractString=".png"
    )
        function callback(; kwargs...)
            mkpath(directory)
            for (name, image) in kwargs
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

    skipnanormissing(arr::AbstractArray) = filter(x -> !ismissing(x) && !isnan(x), arr)
end
