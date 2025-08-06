"""
function segmentation_comparison(
    validated::SegmentedImage, measured::SegmentedImage
)::@NamedTuple{recall::Real, precision::Real, F_score::Real}

Compares two SegmentedImages and returns values describing how similar the segmentations are.

Treats the segment labeled `0` as background.

Measures:
- normalized_{validated,measured}_area: fraction of image covered by segments
- fractional_intersection: fraction of the validated segments covered by measured segments
"""
function segmentation_comparison(
    validated::Union{SegmentedImage}, measured::Union{SegmentedImage}
)::@NamedTuple{recall::Real, precision::Real, F_score::Real}
    validated_binary = binarize(validated)
    measured_binary = binarize(measured)
    intersection = measured_binary .&& validated_binary
    recall = sum(intersection) / sum(validated_binary)
    precision = sum(intersection) / sum(measured_binary)
    F_score = 2 * (precision * recall) / (precision + recall)
    return (; recall, precision, F_score)
end

function binarize(segments::SegmentedImage)
    return labels_map(segments) .> 0
end

function segmentation_comparison(;
    validated::Union{SegmentedImage}, measured::Union{SegmentedImage}
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
