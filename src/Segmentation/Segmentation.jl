module Segmentation

export bwdist

using Images: distance_transform, feature_transform

"""
    bwdist(bwimg)

Distance transform for binary image `bwdist`.
"""
function bwdist(bwimg::AbstractArray{Bool})::AbstractArray{Float64}
    return distance_transform(feature_transform(bwimg))
end

end
