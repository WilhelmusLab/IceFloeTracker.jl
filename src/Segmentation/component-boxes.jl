# TODO: Remove this file once the fixed version of `ImageMorphology.component_boxes` is released.
# Temporary local copy of `ImageMorphology.component_boxes`, carrying the fix in
# JuliaImages/ImageMorphology.jl#146. Delete this file, drop its `include`, and
# restore `component_boxes` to the `import Images` list in regionprops.jl once a
# release containing that fix is available.
#
# Adapted from ImageMorphology.jl (MIT, Copyright (c) 2013-2017: Tim Holy and
# contributors), src/connected.jl.

import OffsetArrays: OffsetArray

"""
    _component_boxes(A)

Minimal bounding box of each label in `A`, as `CartesianIndices`.

`A` must be a labeled array whose minimum is `0` or `1`. When it contains
background the result is shifted to 0-based indexing, so the background region
is entry `0` and label `i` is entry `i`; when its minimum is `1` the result is
indexed `1:maximum(A)` and has no entry `0`. A label absent from `A` gets an
empty box.
"""
function _component_boxes(A::AbstractArray{T,N}) where {T<:Integer,N}
    mn, mx = extrema(A)
    if !(mn == 0 || mn == 1)
        throw(
            ArgumentError(
                "The input labeled array should contain background label `0` as the minimum value",
            ),
        )
    end
    boxes = if mn == 1
        OffsetArray(Matrix{CartesianIndex{N}}(undef, mx, 2), 0, 0)
    elseif mn == 0
        OffsetArray(Matrix{CartesianIndex{N}}(undef, 1 + mx, 2), -1, 0)
    end
    R = CartesianIndices(A)
    boxes[:, 1] .= Ref(last(R))
    boxes[:, 2] .= Ref(first(R))

    @inbounds for i in R
        label_idx = A[i]
        boxes[label_idx, 1] = min(i, boxes[label_idx, 1])
        boxes[label_idx, 2] = max(i, boxes[label_idx, 2])
    end

    return [boxes[i, 1]:boxes[i, 2] for i in axes(boxes, 1)]
end
