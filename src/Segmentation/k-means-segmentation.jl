import Images: Gray, Float64
import StatsBase: StatsBase
import Random: Random
import Clustering: Clustering, assignments

"""
    kmeans_segmentation(gray_image, ice_labels;)

Apply k-means segmentation to a gray image to isolate a cluster group representing sea ice. Returns a binary image with ice segmented from background.

# Arguments
- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

"""
function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}},
    ice_labels::Union{Vector{Int64},BitMatrix,AbstractArray{<:Gray}},
    k::Int64=4,
    maxiter::Int64=50,
)::BitMatrix
    segmented = kmeans_segmentation(gray_image; k=k, maxiter=maxiter)

    ## Isolate ice floes and contrast from background
    segmented_ice = get_segmented_ice(segmented, ice_labels)
    return segmented_ice
end

function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}}; k::Int64=4, maxiter::Int64=50, random_seed::Int64=45
)
    Random.seed!(random_seed)

    ## NOTE(tjd): this clusters into k classes and solves iteratively with a max of maxiter iterations
    feature_classes = Clustering.kmeans(
        vec(gray_image), k; maxiter=maxiter, display=:none, init=:kmpp
    )

    class_assignments = assignments(feature_classes)

    segmented = reshape(class_assignments, size(gray_image))

    return segmented
end

function get_segmented_ice(
    segmented::Matrix{Int64}, ice_labels::Union{Vector{Int64},BitMatrix}
)
    ## Same principle as the get_nlabels function
    ## Has the weakness that only one segment can ever be chosen.
    isempty(ice_labels) && return falses(size(segmented))
    return segmented .== StatsBase.mode(segmented[ice_labels])
end

function get_segmented_ice(segmented::Matrix{Int64}, ice_labels::AbstractArray{<:Gray})
    boolean_map = ice_labels .|> Bool
    !(any(boolean_map)) && return falses(size(segmented))
    return segmented .== StatsBase.mode(segmented[boolean_map])
end
