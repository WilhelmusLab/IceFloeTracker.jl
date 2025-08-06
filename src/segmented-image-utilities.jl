"""
segmentation_comparison(
    validated::Union{SegmentedImage,Nothing},
    measured::Union{SegmentedImage,Nothing},
)::NamedTuple

segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}=nothing,
    measured::Union{SegmentedImage,Nothing}=nothing,
)::NamedTuple

Compares two SegmentedImages and returns values describing how similar the segmentations are.

Treats the segment labeled `0` as background.

Measures:
- normalized_{validated,measured}_area: fraction of image covered by segments
- fractional_intersection: fraction of the validated segments covered by measured segments
"""
function segmentation_comparison(
    validated::Union{SegmentedImage,Nothing}, measured::Union{SegmentedImage,Nothing}
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

function segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}=nothing,
    measured::Union{SegmentedImage,Nothing}=nothing,
)::NamedTuple
    return segmentation_comparison(validated, measured)
end
