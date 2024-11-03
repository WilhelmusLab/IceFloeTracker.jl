function imcomplement(img::Matrix{T}) where {T<:Union{Unsigned,Int}}
    return 255 .- img
end

function imcomplement(img::Matrix{Gray{Float64}})
    return 1 .- img
end

function reconstruct(img, se, type, invert::Bool=true)
    if type == "dilation"
        morphed = IceFloeTracker.MorphSE.dilate(img, se)
    elseif type == "erosion"
        morphed = IceFloeTracker.MorphSE.erode(img, se)
    else
        throw(ArgumentError("Invalid type: $type. Must be 'dilation' or 'erosion'."))
    end

    if invert
        morphed = imcomplement(morphed)
        img = imcomplement(img)
    end

    return IceFloeTracker.MorphSE.mreconstruct(IceFloeTracker.MorphSE.dilate, morphed, img)
end

function reconstruct_dilation(img, se)
    return reconstruct(img, se, "dilation", true)
end

function reconstruct_erosion(img, se)
    return reconstruct(img, se, "erosion", false)
end
