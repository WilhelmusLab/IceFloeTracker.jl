# Helper functions

"""
    bwdist(bwimg)

Distance transform for binary image `bwdist`.
"""
function bwdist(bwimg::AbstractArray{Bool})::AbstractArray{Float64}
    return Images.distance_transform(Images.feature_transform(bwimg))
end
