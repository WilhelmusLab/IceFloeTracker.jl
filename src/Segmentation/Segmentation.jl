module Segmentation

export IceFloeSegmentationAlgorithm,
    find_ice_labels,
    get_ice_labels_mask,
    get_ice_labels,
    get_ice_masks,
    find_ice_mask,
    tiled_adaptive_binarization,
    IceDetectionAlgorithm,
    IceDetectionThresholdMODIS721,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionLopezAcosta2019,
    regionprops_table,
    regionprops,
    SegmentationComparison,
    segmentation_comparison,
    SegmentationSummary,
    segmentation_summary,
    binarize_segments,
    kmeans_segmentation,
    addlatlon!,
    convertcentroid!,
    converttounits!,
    get_ice_peaks

include("abstract-algorithms.jl")
include("find-ice-labels.jl")
include("ice-masks.jl")
include("k-means-segmentation.jl")
include("regionprops.jl")
include("segmented-image-utilities.jl")

end
