module Segmentation

export find_ice_labels,
    get_ice_labels_mask,
    get_ice_masks,
    tiled_adaptive_binarization,
    IceDetectionAlgorithm,
    IceDetectionThresholdMODIS721,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionLopezAcosta2019,
    regionprops_table,
    SegmentationComparison,
    segmentation_comparison,
    SegmentationSummary,
    segmentation_summary,
    binarize_segments

include("find_ice_labels.jl")
include("ice_masks.jl")
include("regionprops.jl")
include("segmented-image-utilities.jl")

end
