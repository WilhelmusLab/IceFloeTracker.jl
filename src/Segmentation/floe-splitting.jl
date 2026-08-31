"""
Functions for separating ice floes in binary images.
"""

import Images: 
    distance_transform,
    feature_transform,
    label_components,
    opening,
    component_indices,
    imfill
import ..Morphology:
    strel_disk

"""
    dist_morph_split(
        binary_floes::BitMatrix;
        min_floe_size::Int64=64,
        max_hole_fill::Int64=2000,
        max_depth::Int64=5,
        max_expand::Int64=3,
        strel=strel_disk(3)
    )

Method to split objects in a binary image using image morphology and the distance transform. The algorithm
operates by calculating the distance transform, which computes the distance from each labeled pixel to the background.
There are two steps: creating a ``pyramid'', then stepping down from the top of the pyramid and re-labeling or expanding
shapes as needed.

For each distance d up to `max_distance`, select pixels that are greater than that distance. Perform morphological opening,
fill holes up to `max_hole_fill`, then label components. Each of these layers is a level in the pyramid.

Then, starting from the highest level of the pyramid, check to see whether objects in the next layer down contain multiple
objects in the current layer. If an object at layer ``d-1`` contains only object at layer ``d``, then keep the object at layer ``d-1``.
Otherwise, expand the labels by `max_expand`, then intersect the expanded labels with the containing object at layer ``d-1``.

After traversing the pyramid, relabel matrix, and remove any objects smaller than the `min_floe_size`.

"""
function dist_morph_split(
    binary_floes::BitMatrix;
    max_hole_fill::Int64=2000,
    max_depth::Int64=5,
    max_depth_ratio::Real=0.3,
    max_expand::Int64=3,
    opening_strel=strel_disk(3),
)
    dist = distance_transform(feature_transform(.!binary_floes))
    # Initialize with one run of opening
    levels = Dict(0 => label_components(opening(dist .> 0, opening_strel)))

    ### Build pyramid - each size is the opened and filled thresholded image for a given distance
    for dist_threshold in 1:max_depth
        markers = opening(dist .> dist_threshold, opening_strel)
        markers .= .!imfill(.!markers, (0, max_hole_fill))
        labeled_markers = label_components(markers)
        maximum(labeled_markers) == 0 && break

        labels = filter(r -> r != 0, unique(labeled_markers))
        indices = component_indices(labeled_markers)
        
        # check 1: Remove components with no intersection with the layer below
        remove_list = _nonoverlapping_labels(levels[dist_threshold - 1], indices, labels)
        _remove_labels!(labeled_markers, indices, remove_list)
        filter!(r -> r ∉ remove_list, labels)

        # check 2: Remove components which fail the max_depth_ratio to component maximum depth test
        maximum_depths = Dict(L => maximum(dist[indices[L]]) for L in labels)
        remove_list = [L for L ∈ labels if max_depth_ratio * maximum_depths[L] < dist_threshold]
        _remove_labels!(labeled_markers, indices, remove_list)
        levels[dist_threshold] = labeled_markers
    end
    max_depth = maximum([d for d in keys(levels)])
    final_labels = copy(levels[max_depth])

    ### Descend pyramid
    for dist_threshold in max_depth:-1:1
        # Get indices from level d-1
        indices = component_indices(levels[dist_threshold - 1])
        labels = filter(r -> r != 0, unique(levels[dist_threshold - 1]))
        # Expand indices at level d
        expanded = expand_labels(levels[dist_threshold], max_expand)
        for L in labels
            matched_labels = unique(levels[dist_threshold][indices[L]])
            # If intersection of the label at level
            if (0 ∈ matched_labels) && (length(matched_labels) <= 2)
                final_labels[indices[L]] .= L
                continue
            end
            # Otherwise, expand the current level, and set the next level down to the expanded indices.
            # May need to check the number of matched labels in the expanded image.
            levels[dist_threshold - 1][indices[L]] .= expanded[indices[L]]
            final_labels[indices[L]] .= expanded[indices[L]]
        end
    end
    return label_components(final_labels)
end

function _nonoverlapping_labels(other, indices, labels)
    return [
        label for label in labels
        if maximum(other[indices[label]]) == 0
    ]
end

function _remove_labels!(output, indices, remove_labels)
    for L in remove_labels
        if L != 0
            output[indices[L]] .= 0
        end
    end
end
