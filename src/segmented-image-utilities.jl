"""
Results of a segmentation comparison
"""
SegmentationComparison = @NamedTuple begin
    recall::Union{Real,Missing}
    precision::Union{Real,Missing}
    F_score::Union{Real,Missing}
end

"""
function segmentation_comparison(
    validated::SegmentedImage, measured::SegmentedImage
)::@NamedTuple{recall::Real, precision::Real, F_score::Real}

Compares two SegmentedImages and returns values describing how similar the segmentations are.

This treats the segment labeled `0` as background.

Measures:
- precision: rate at which pixels in `validated` segments belong to `measured` segments
- recall: rate at which pixels in `measured` segments belong to `validated` segments
- F_score: harmonic mean of precision and recall
"""
function segmentation_comparison(
    validated::Union{SegmentedImage,Nothing}, measured::Union{SegmentedImage,Nothing}
)::SegmentationComparison
    (isnothing(validated) || isnothing(measured)) &&
        return (; recall=missing, precision=missing, F_score=missing)

    validated_binary = binarize(validated)
    measured_binary = binarize(measured)
    intersection = measured_binary .&& validated_binary
    recall = sum(intersection) / sum(validated_binary)
    precision = sum(intersection) / sum(measured_binary)
    F_score = 2 * (precision * recall) / (precision + recall)
    return (; recall, precision, F_score)
end

function segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}, measured::Union{SegmentedImage,Nothing}
)::SegmentationComparison
    return segmentation_comparison(validated, measured)
end

function binarize(segments::SegmentedImage)
    return labels_map(segments) .> 0
end

"""
Results of a segmentation comparison
"""
SegmentationSummary = @NamedTuple begin
    labeled_fraction::Real
end

function segmentation_summary(segmented::SegmentedImage)::SegmentationSummary
    binary = binarize(segmented)
    non_zero_area = sum(binary)
    labeled_fraction = non_zero_area / length(binary)
    return (; labeled_fraction)
end
