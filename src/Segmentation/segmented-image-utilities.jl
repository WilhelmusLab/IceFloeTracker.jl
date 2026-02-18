import ImageSegmentation: SegmentedImage, labels_map, segment_labels, segment_mean, segment_pixel_count
using DataFrames 
using Random

"""
Results of a segmentation comparison
"""
SegmentationComparison = @NamedTuple begin
    recall::Union{Real,Missing}
    precision::Union{Real,Missing}
    F_score::Union{Real,Missing}
end

"""
    segmentation_comparison(
        validated::SegmentedImage, measured::SegmentedImage
    )::@NamedTuple{recall::Real, precision::Real, F_score::Real}

Compares two SegmentedImages and returns values describing how similar the segmentations are.

This treats the segment labeled `0` as background.

Measures:
- precision: rate at which pixels in `validated` segments belong to `measured` segments
- recall: rate at which pixels in `measured` segments belong to `validated` segments
- F_score: harmonic mean of precision and recall
"""
function segmentation_comparison(
    validated::Union{SegmentedImage,Nothing}, measured::Union{SegmentedImage,Nothing}
)::SegmentationComparison
    (isnothing(validated) || isnothing(measured)) &&
        return (; recall=missing, precision=missing, F_score=missing)

    validated_binary = binarize_segments(validated)
    measured_binary = binarize_segments(measured)
    intersection = @. Bool(measured_binary) && Bool(validated_binary)
    recall = sum(intersection) / sum(validated_binary)
    precision = sum(intersection) / sum(measured_binary)
    F_score = 2 * (precision * recall) / (precision + recall)
    return (; recall, precision, F_score)
end

function segmentation_comparison(;
    validated::Union{SegmentedImage,Nothing}, measured::Union{SegmentedImage,Nothing}
)::SegmentationComparison
    return segmentation_comparison(validated, measured)
end

"""
Results of a segmentation comparison
"""
SegmentationSummary = @NamedTuple begin
    labeled_fraction::Union{Real,Missing}
    segment_count::Union{Real,Missing}
end

function segmentation_summary(segmented::Union{SegmentedImage,Nothing})::SegmentationSummary
    if !isnothing(segmented)
        binary = binarize_segments(segmented)
        non_zero_area = sum(binary)
        labeled_fraction = non_zero_area / length(binary)
        segment_count = length(segment_labels(segmented))
        return (; labeled_fraction, segment_count)
    else
        return (; labeled_fraction=missing, segment_count=missing)
    end
end

"""
    binarize_segments(segments::SegmentedImage)

Find pixels in a segmented image with non-zero labels and return as a grayscale image.
"""
function binarize_segments(segments::SegmentedImage)::AbstractArray{Gray}
    return Gray.(labels_map(segments) .> 0)
end


"""
    stitch_clusters(tiles, segmented_image, minimum_overlap, grayscale_threshold)

Stitches clusters across tile boundaries based on neighbor with largest shared boundary.
The algorithm finds all pairs of segment labels at the tile edges. Then, we count the number of 
times each right-hand label is paired to a left-hand label, and for pairs with at least `minimum_overlap` pixel overlap,
the right-hand label is assigned as a candidate pair to the left-hand label. If the difference in grayscale
intensity is less than `grayscale_threshold`, the objects are merged. The function returns an image index map.
"""
function stitch_clusters(segmented_image, tiles, minimum_overlap=4, grayscale_threshold=0.1) 
    grayscale_magnitude(c) = Float64(Gray(c))
    
    idxmap = deepcopy(labels_map(segmented_image))
    n, m = size(idxmap)
    
    for tile in tiles
        nrange, mrange = tile
        tn = maximum(nrange)
        tm = maximum(mrange)
        label_pairs = []
        if tn != n
            push!(label_pairs, vec([(x, y) for (x, y) in zip(idxmap[tn, :], idxmap[tn .+ 1, :])]))
        end
        
        if tm != m
            push!(label_pairs, vec([(x, y) for (x, y) in zip(idxmap[:, tm], idxmap[:, tm .+ 1])]))
        end

        if !isempty(label_pairs)

            # create a dataframe with the results
            label_pairs = vcat(label_pairs...)
            label_pairs = reshape(reinterpret(Int64, label_pairs), (2,:))
            df = DataFrame(left=label_pairs[1,:], right=label_pairs[2,:])
            
            # groupby right -> left pairs and get counts
            df_counts = combine(groupby(df, [:right, :left]), nrow => :count)
            
            # only use pairs that overlap by at least 2 pixels.
            df_counts = df_counts[df_counts.count .>= minimum_overlap, :]

            # only use pairs where left is not equal to right
            df_counts = df_counts[df_counts.right .!= df_counts.left, :]

            # don't merge if the segments are too different in color
            left_brightness = [grayscale_magnitude(segment_mean(segmented_image, l)) for l in df_counts.left]
            right_brightness = [grayscale_magnitude(segment_mean(segmented_image, r)) for r in df_counts.right]
            diff_means = abs.(right_brightness .- left_brightness)
            df_counts = df_counts[diff_means .< grayscale_threshold, :]
            
            if !isempty(df_counts)
                # now find the maximum overlapping segment for each
                df_pairs = combine(sdf -> sdf[argmax(sdf.count), [:right, :left, :count]], groupby(df_counts, :right))
            
                # make a lookup table and lookup function
                lut = Dict(ri => li for (ri, li) in zip(df_pairs.right, df_pairs.left))
                
                lookup_remap(ii) = begin
                    if haskey(lut, ii)
                       return lut[ii]
                    end
                    return ii
                end
    
                idxmap .= map(i -> lookup_remap(i), idxmap)
            end
        end    
    end
    return idxmap
end

"""
    _get_random_color(seed)

Convenience function to produce a random RGB color.
"""
function _get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

"""
    view_seg_random(s::SegmentedImage)

Produce an array with the segment mean mapped to each segment label.
If the SegmentedImage was produced with color type (e.g., RGB, Gray), then
the result will be an image.
"""
function view_seg(s)
    map(i->segment_mean(s,i), labels_map(s))
end

"""
    view_seg_random(s::SegmentedImage)

Produce an RGB image with a random color for each unique segment in `s`.
"""
function view_seg_random(s)
    map(i->_get_random_color(i), labels_map(s))
end
