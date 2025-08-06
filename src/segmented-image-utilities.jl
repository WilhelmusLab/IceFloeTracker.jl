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
        normalized_validated_area = segmentation_summary(validated).normalized_non_zero_area
    else
        normalized_validated_area = missing
    end

    if !isnothing(measured)
        normalized_measured_area = segmentation_summary(measured).normalized_non_zero_area
    else
        normalized_measured_area = missing
    end

    if !isnothing(validated) && !isnothing(measured)
        intersection = Gray.(channelview(measured_binary) .&& channelview(validated_binary))
        recall = sum(channelview(intersection)) / validated_area
        precision = sum(channelview(intersection)) / measured_area
        F_score = 2 * (precision * recall) / (precision + recall)
    else
        recall = missing
        precision = missing
    end

    return (;
        normalized_validated_area,
        normalized_measured_area,
        fractional_intersection=recall,
        recall,
        precision,
        F_score,
    )
end

function segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}=nothing,
    measured::Union{SegmentedImage,Nothing}=nothing,
)::NamedTuple
    return segmentation_comparison(validated, measured)
end

function segmentation_summary(
    image::SegmentedImage
)::@NamedTuple{normalized_non_zero_area::Real}
    image_binary = Gray.(labels_map(image) .> 0)
    non_zero_area = sum(channelview(image_binary))
    normalized_non_zero_area = non_zero_area / length(channelview(image_binary))
    return (; normalized_non_zero_area)
end
