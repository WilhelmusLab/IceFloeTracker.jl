import Images: Gray, Float64, SegmentedImage
import StatsBase: StatsBase
import Random: Random
import Clustering: Clustering, assignments, kmeans
import ..ImageUtils: @tilewise, tile_index_map

"""
    kmeans_segmentation(gray_image; k=4, maxiter=50, random_seed=45)

Wrapper for Clustering.kmeans which accepts a grayscale image and returns a SegmentedImage object. Optionally,
one can specify the number of clusters `k``, the maximum number of iterations `maxiter`, and the seed for the
random number generator, `random_seed`. Returns a SegmentedImage object.
"""
function kmeans_segmentation(
    gray_image::AbstractArray{<:AbstractGray};
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
)
    Random.seed!(random_seed)

    feature_classes = kmeans(vec(gray_image), k; maxiter=maxiter, display=:none, init=:kmpp)

    class_assignments = assignments(feature_classes)

    segmented = reshape(class_assignments, size(gray_image))

    return SegmentedImage(gray_image, segmented)
end

function kmeans_segmentation(
    gray_image::AbstractArray{<:AbstractGray},
    tiles;
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
    minimum_overlap=4,
    grayscale_threshold=0.1,
)
    # Each tile is clustered on its own, so every tile has labels 1:k. Shifting the labels
    # of the i-th tile by k*i makes them unique across tiles: tile 1 gets k+1:2k, tile 2
    # gets 2k+1:3k, and so on. Pixels outside every tile stay 0.
    segment(img) = _kmeans_labels(img; k=k, maxiter=maxiter, random_seed=random_seed)
    canvas = zeros(Int64, size(gray_image))
    offsets = k .* tile_index_map(tiles, size(gray_image))

    indexmap = @tilewise tiles segment($gray_image) into = canvas # labels 1:k per tile
    indexmap .+= offsets # unique labels across tiles

    # Merge clusters that continue across tile edges
    indexmap .= labels_map(
        stitch_clusters(
            SegmentedImage(gray_image, indexmap),
            tiles,
            minimum_overlap,
            grayscale_threshold,
        ),
    )

    return SegmentedImage(gray_image, indexmap)
end

"""
    _kmeans_labels(gray_image; k, maxiter, random_seed)

Label map of the k-means segmentation of `gray_image`, with labels `1:k`.
"""
_kmeans_labels(gray_image; k, maxiter, random_seed) =
    labels_map(kmeans_segmentation(gray_image; k=k, maxiter=maxiter, random_seed=random_seed))

"""
    kmeans_binarization(gray_image, false_color_image; kwargs...)

Produce a binarized image by identifying pixels of bright ice, performing k-means clustering, and then selecting the k-means cluster
containing the largest fraction of bright ice pixels. If no bright ice pixels are detected, then a blank matrix is returned.

# Positional Arguments
- `gray_image`: Grayscale image to segment using k-means.
- `falsecolor_image`: MODIS 7-2-1 falsecolor image, to be sent to the specified `ice_labels_algorithm`. It is recommended that this image be landmasked.

# Keyword arguments
- `ice_labels_algorithm`: Binarization function to find sea ice pixels
- `k`: Number of k-means clusters
- `maxiter`: Maximum number of iterations for k-means algorithm
- `random_seed`: Seed for the random number generator
"""
function kmeans_binarization(
    gray_image,
    falsecolor_image;
    cluster_selection_algorithm=IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01),
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
    threshold=1, # TODO: Make the test FirstNonZero algo more robust, case 14 succeeds with only 1 or 2 pixels which is not stable.
)::BitMatrix
    selected_labels = cluster_selection_algorithm(falsecolor_image) .> 0
    isempty(selected_labels) && return falses(size(gray_image))
    sum(selected_labels) < threshold && return falses(size(gray_image))

    segmented = kmeans_segmentation(
        gray_image; k=k, maxiter=maxiter, random_seed=random_seed
    )

    return segmented.image_indexmap .==
        StatsBase.mode(segmented.image_indexmap[selected_labels])
end

"""
    kmeans_binarization(gray_image, false_color_image, tiles; kwargs...)

Produce a binarized image tilewise by identifying pixels of bright ice, performing k-means clustering, and then selecting the k-means cluster
containing the largest fraction of bright ice pixels in each tile. If less than a threshold number of bright ice pixels are detected, then a blank matrix is returned.


# Positional Arguments
- `gray_image`: Grayscale image to segment using k-means.
- `falsecolor_image`: MODIS 721 false color image to be supplied to the cluster selection algorithm.
- `tiles`: Tiled iterator (e.g. from `IceFloeTracker.get_tiles()`)

# Keyword arguments
- `k`: Number of k-means clusters. Default 4.
- `maxiter`: Maximum number of iterations for k-means algorithm. Default 50.
- `random_seed`: Seed for the random number generator. Default 45.
- `cluster_selection_algorithm`: Binarization function to find the k-means cluster to set to 1; all other clusters set to 0.
- `threshold`: Minimum number of ice pixels to trigger selection of k cluster
- `minimum_overlap`: Argument to `stitch_clusters`, minimum number of pixels on boundary for merge
- `grayscale_threshold`: Argument to `stitch_clusters`, maximum grayscale difference for merge
"""
function kmeans_binarization(
    gray_image,
    falsecolor_image,
    tiles;
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
    cluster_selection_algorithm,
    ice_pixels_threshold=1,
    minimum_overlap=4,
    grayscale_threshold=0.1,
)::BitMatrix
    kseg = kmeans_segmentation(
        gray_image,
        tiles;
        k=k,
        random_seed=random_seed,
        maxiter=maxiter,
        minimum_overlap=minimum_overlap,
        grayscale_threshold=grayscale_threshold,
    )
    labels = labels_map(kseg)
    ice = cluster_selection_algorithm(falsecolor_image) .> 0

    select_cluster(labels_tile, ice_tile) = _select_ice_cluster(
        labels_tile, ice_tile; ice_pixels_threshold=ice_pixels_threshold
    )
    canvas = falses(size(gray_image))

    return @tilewise tiles select_cluster($labels, $ice) into = canvas
end

"""
    _select_ice_cluster(labels, ice; ice_pixels_threshold)

Return a `BitMatrix` that is true where `labels` equals the most common label among the `ice`
pixels. If `ice` has `ice_pixels_threshold` or fewer true pixels, return all false.
"""
function _select_ice_cluster(labels, ice; ice_pixels_threshold)
    sum(ice) > ice_pixels_threshold || return falses(size(ice))
    return labels .== StatsBase.mode(labels[ice])
end
