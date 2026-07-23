module Segmentation

export addlatlon!,
    BenkridCrookes,
    binarize_segments,
    component_convex_areas,
    component_floes,
    component_perimeters,
    convertcentroid!,
    converttounits!,
    converttounits,
    ConvexAreaEstimationAlgorithm,
    expand_labels,
    get_ice_peaks,
    kmeans_binarization,
    kmeans_segmentation,
    tiled_adaptive_binarization,
    IceDetectionAlgorithm,
    IceDetectionThresholdMODIS721,
    IceDetectionBrightnessMidpoint,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    PerimeterEstimationAlgorithm,
    PolygonConvexArea,
    PixelConvexArea,
    regionprops_table,
    regionprops,
    remove_small_segments!,
    remove_large_segments!,
    remove_low_contrast_segments!,
    segment_mean_map,
    SegmentationComparison,
    segmentation_comparison,
    SegmentationSummary,
    segmentation_summary,
    stitch_clusters,
    view_seg,
    view_seg_random

include("ice-detection-algorithms.jl")
include("k-means-segmentation.jl")
include("regionprops.jl")
include("segmented-image-utilities.jl")

end
