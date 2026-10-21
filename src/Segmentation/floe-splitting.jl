"""
Functions for separating ice floes in binary images.
"""

import Images: 
    distance_transform,
    feature_transform,
    label_components,
    opening,
    component_indices,
    imfill,
    watershed

import ..Morphology:
    strel_disk,
    fill_holes!

"""
    dist_morph_split(
        binary_floes::BitMatrix;
        min_floe_size::Int64=64,
        max_hole_fill::Int64=2000,
        max_depth::Int64=5,
        max_expand::Int64=2,
        strel=strel_disk(3)
    )

Method to split objects in a binary image using image morphology and the distance transform. The algorithm
operates by calculating the distance transform, which computes the (minimum) distance from each labeled pixel to the background.
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
    min_floe_size::Int64=64,
    max_hole_fill::Int64=2000,
    max_depth::Int64=25,
    max_depth_ratio::Real=0.7,
    max_expand::Int64=3,
    opening_strel=strel_disk(3),
)::Matrix{Int64}
    dist = distance_transform(feature_transform(.!binary_floes))
    # Initialize level 0, exit early if opening removes all potential floes
    markers = opening(dist .> 0, opening_strel)
    markers .= .!imfill(.!markers, (0, max_hole_fill))
    labeled_markers = label_components(markers)
    maximum(labeled_markers) == 0 && return markers 
    levels = Dict(0 => labeled_markers)
   
    ### Build pyramid - each size is the opened and filled thresholded image for a given distance
    for dist_threshold in 1:max_depth
        markers = opening(dist .> dist_threshold, opening_strel)
        markers .= .!imfill(.!markers, (0, max_hole_fill))
        labeled_markers = label_components(markers)
        # maximum(labeled_markers) == 0 && break

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
    remove_small_segments!(final_labels, min_floe_size)
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

# ebseg python code
import IceFloeTracker.Morphology: fill_holes

"""

EE method from Buckley et al. 2024

"""
function erosion_expansion(
        gray_img, # image to use in the watershed transformation
        binary_img; # bitmatrix where ice is marked true
        itmin=1, # minimum number of iterations
        itmax=8, # maximum number of iterations
        step=1, # spacing between iterations
        strel=strel_diamond((3,3)),
        minimum_area=16,
        save_figs=true
        )

    input = copy(binary_img)
    input_no = copy(binary_img)
    
    output = zeros(Int64, size(binary_img))
    label_offset = 0
    
    for iters in itmax:-step:itmin
        eroded_ice_mask = copy(binary_img)

        # Note that copy is needed for in-place morphological operations, since
        # the in-place operation will delete the source image along the way.
        
        foreach(_ -> erode!(eroded_ice_mask, copy(eroded_ice_mask), strel), 1:iters)
        fill_holes!(eroded_ice_mask, copy(eroded_ice_mask), strel_diamond((3,3)))
        markers = label_components(eroded_ice_mask) .+ 1
        dilated_ice_mask = copy(eroded_ice_mask)
        foreach(_ -> dilate!(dilated_ice_mask, copy(dilated_ice_mask), strel), 1:iters)

        # Select the pixels between the dilated and the eroded sections.
        # These are where we'll re-grow floes.
        # The phrasing "unknown" is based on the CV2 watershed description
        # (e.g. https://docs.opencv.org/4.13.0/d3/db4/tutorial_py_watershed.html)
        # Setting the background to 1 helps keep it separate at the end.
        
        unknown = copy(dilated_ice_mask)
        unknown[eroded_ice_mask .> 0] .= 0
        markers[unknown .> 0] .= 0
        
        # Dilating the markers replaces values with the maximum in the neighborhood,
        # which can have surprising results with randomly-labeled data.
        foreach(_ -> dilate!(markers, copy(markers), strel), 1:iters)
            
        # Use the watershed transformation to recover some detail of floe edges
        watershed_labels = labels_map(watershed(complement.(gray_img), markers))
        watershed_labels .-= 1
        
        # Set anything originally classified as water to background
        watershed_labels[binary_img .== 0] .= 0
        
        # Clear objects intersecting the land mask (TBD - no land in this image)

        # Clear objects where we've already got output labels
        watershed_labels[output .> 0] .= 0

        # Re-label since regions often are split up at this stage
        watershed_labels .= label_components(watershed_labels)
        
        # Remove small segments
        remove_small_segments!(watershed_labels, maximum([minimum_area, iters^4]))

        if save_figs
            fname = "identification_round_"*string(iters)*".png"
            save(fname, view_seg_random(SegmentedImage(binary_img, watershed_labels)))
        end

        # Input for next stage:
        # backround in the watershed, positive in the ice mask, positive at beginning of round
        input .= (watershed_labels .== 0) .&& (input .> 0) .&& (binary_img .> 0)

        # Update the output
        new_label_mask = watershed_labels .> 0
        output[new_label_mask] .= watershed_labels[new_label_mask] .+ label_offset
        label_offset = maximum(output)
    end
    output .= label_components(output)
    
    # add min brightness step (TBD)
    
    return output
end
