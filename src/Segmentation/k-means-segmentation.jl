import Images: Gray, Float64, SegmentedImage
import StatsBase: StatsBase
import Random: Random
import Clustering: Clustering, assignments, kmeans

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
    random_seed::Int64=45
)
    Random.seed!(random_seed)

    feature_classes = kmeans(
        vec(gray_image), k; maxiter=maxiter, display=:none, init=:kmpp
    )

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
    grayscale_threshold=0.1
)

    indexmap = zeros(Int64, size(gray_image))
    for (offset, tile) in enumerate(tiles)
        indexmap[tile...] .= labels_map(kmeans_segmentation(gray_image[tile...]; k=k, maxiter=maxiter, random_seed=random_seed))
        indexmap[tile...] .+= k * offset
        offset += 1
    end

    indexmap .= stitch_clusters(SegmentedImage(gray_image, indexmap), tiles, minimum_overlap, grayscale_threshold) 
    return SegmentedImage(gray_image, indexmap)
end


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
    cluster_selection_algorithm=IceDetectionBrightnessPeaksMODIS721(0.2, 0.3),
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
    threshold=1 # TODO: Make the test FirstNonZero algo more robust, case 14 succeeds with only 1 or 2 pixels which is not stable.
)::BitMatrix

    selected_labels = cluster_selection_algorithm(falsecolor_image) .> 0
    isempty(selected_labels) && return falses(size(gray_image))
    sum(selected_labels) < threshold && return falses(size(gray_image))

    segmented = kmeans_segmentation(gray_image; k=k, maxiter=maxiter, random_seed=random_seed)

    return segmented.image_indexmap .== StatsBase.mode(segmented.image_indexmap[selected_labels])
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
    grayscale_threshold=0.1
)::BitMatrix
    out = falses(size(gray_image))
    kseg = kmeans_segmentation(gray_image, tiles;
                k=k, random_seed=random_seed,
                maxiter=maxiter,
                minimum_overlap=minimum_overlap,
                grayscale_threshold=grayscale_threshold)
    labels = labels_map(kseg)
    ice = cluster_selection_algorithm(falsecolor_image) .> 0
        
    for tile in tiles
        if sum(ice[tile...]) > ice_pixels_threshold
            k = StatsBase.mode(labels[tile...][ice[tile...]])
            out[tile...] .= labels[tile...] .== k
        end
    end

    return out
end