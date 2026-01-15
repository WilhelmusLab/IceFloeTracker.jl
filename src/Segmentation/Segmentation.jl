module Segmentation

export
    addlatlon!,
    BenkridCrookes,
    binarize_segments,
    component_floes,
    component_perimeters,
    convertcentroid!,
    converttounits!,
    ConvexAreaEstimationAlgorithm,
    find_ice_labels,
    get_ice_labels_mask,
    get_ice_labels,
    get_ice_masks,
    get_ice_peaks,
    find_ice_mask,
    kmeans_binarization,
    kmeans_segmentation,
    tiled_adaptive_binarization,
    IceDetectionAlgorithm,
    IceDetectionThresholdMODIS721,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionLopezAcosta2019,
    IceFloeSegmentationAlgorithm,
    PerimeterEstimationAlgorithm,
    PolygonConvexArea,
    PixelConvexArea,
    regionprops_table,
    regionprops,
    SegmentationComparison,
    segmentation_comparison,
    SegmentationSummary,
    segmentation_summary,
    stitch_clusters

include("abstract-algorithms.jl")
include("find-ice-labels.jl")
include("ice-masks.jl")
include("k-means-segmentation.jl")
include("regionprops.jl")
include("segmented-image-utilities.jl")

end
