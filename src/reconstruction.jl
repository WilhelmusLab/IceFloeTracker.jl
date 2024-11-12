function reconstruct(img, se, type, invert::Bool=true)
    !(type == "dilation" || type == "erosion") &&
        throw(ArgumentError("Invalid type: $type. Must be 'dilation' or 'erosion'."))

    type == "dilation" && (morphed = IceFloeTracker.MorphSE.dilate(img, se))
    type == "erosion" &&
        (morphed = IceFloeTracker.sk_morphology.erosion(img; footprint=collect(se)))

    invert && (morphed = imcomplement(morphed); img = imcomplement(img))

    type == "dilation" && return IceFloeTracker.MorphSE.mreconstruct(
        IceFloeTracker.MorphSE.dilate, morphed, img
    )

    return IceFloeTracker.sk_morphology.reconstruction(morphed, img)
end

function reconstruct_erosion(img, se)
    return reconstruct(img, se, "dilation")
end
