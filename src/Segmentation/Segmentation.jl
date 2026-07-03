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
    find_ice_labels,
    find_ice_mask,
    get_ice_labels,
    get_ice_peaks,
    kmeans_binarization,
    kmeans_segmentation,
    tiled_adaptive_binarization,
    IceDetectionAlgorithm,
    IceDetectionThresholdMODIS721,
    IceDetectionBrightnessPeaksMODIS134,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
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
    stitch_clusters,
    view_seg,
    view_seg_random,
    segment_mean_map

include("abstract-algorithms.jl")
include("find-ice-labels.jl")
include("k-means-segmentation.jl")
include("regionprops.jl")
include("segmented-image-utilities.jl")

end
