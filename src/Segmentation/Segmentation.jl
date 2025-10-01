module Segmentation

export bwdist

"""
    dummy_segmentation_function()

Example function for segmentation module

!!! todo "Delete this once real functions are added"
    This function should be removed when real segmentation functions are moved into this module.

"""
function dummy_segmentation_function()
    return "This is a dummy segmentation function."
end

"""
    bwdist(bwimg)

Distance transform for binary image `bwdist`.
"""
function bwdist(bwimg::AbstractArray{Bool})::AbstractArray{Float64}
    return distance_transform(feature_transform(bwimg))
end

end
