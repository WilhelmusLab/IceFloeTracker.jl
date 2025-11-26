import Images: Gray, Float64, SegmentedImage
import StatsBase: StatsBase
import Random: Random
import Clustering: Clustering, assignments, kmeans

"""
    kmeans_segmentation(img; k=4, maxiter=50, random_seed=45, k_offset=0)

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
    ice_labels_algorithm,
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
    ice_labels_threshold=5
)::BitMatrix

    ice_labels = ice_labels_algorithm(falsecolor_image) .> 0
    isempty(ice_labels) && return falses(size(gray_image))
    sum(ice_labels) < ice_labels_threshold && return falses(size(gray_image))

    segmented = kmeans_segmentation(gray_image; k=k, maxiter=maxiter, random_seed=random_seed)

    return segmented.image_indexmap .== StatsBase.mode(segmented.image_indexmap[ice_labels])
end

"""
    kmeans_binarization(gray_image, false_color_image, tiles; kwargs...)

Produce a binarized image tilewise by identifying pixels of bright ice, performing k-means clustering, and then selecting the k-means cluster
containing the largest fraction of bright ice pixels. If no bright ice pixels are detected, then a blank matrix is returned.

Warning: Tilewise processing may result in discontinuities at tile boundaries.

# Positional Arguments
- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

# Keyword arguments
- `k`: Number of k-means clusters
- `maxiter`: Maximum number of iterations for k-means algorithm
- `random_seed`: Seed for the random number generator
- `ice_labels_algorithm`: Binarization function to find sea ice pixels
"""
function kmeans_binarization(
    gray_image,
    falsecolor_image,
    tiles;
    k::Int64=4,
    maxiter::Int64=50,
    random_seed::Int64=45,
    ice_labels_algorithm
)::BitMatrix

    out = falses(size(gray_image))

    # This part will be updated soon. This method results in discontinuities at tile boundaries.
    # To avoid this, we can stitch segmented images first, so that boundary clusters are relabled
    # to match, and then the binarization can be applied using the mode for each tile, allowing overlaps.
    
    for tile in tiles
        out[tile...] .= kmeans_binarization(gray_image[tile...], falsecolor_image[tile...];
                            k=k,
                            maxiter=maxiter,
                            random_seed=random_seed,
                            ice_labels_algorithm=ice_labels_algorithm)
    end
    return out
end
