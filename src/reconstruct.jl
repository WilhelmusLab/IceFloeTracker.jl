# WIP
function imcomplement(img::Matrix{T}) where T<:Union{Unsigned, Int}
    return 255 .- img
end

function imcomplement(img::Matrix{Gray{Float64}})
    return 1 .- img
end

function open_by_reconstruction(img, se)
    marker = IceFloeTracker.MorphSE.erode(img, se)
    return IceFloeTracker.MorphSE.mreconstruct(IceFloeTracker.MorphSE.dilate, marker, img)
end
