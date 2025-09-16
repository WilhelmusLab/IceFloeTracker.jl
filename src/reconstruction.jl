"""
    reconstruct(img, se, type, invert)

Perform closing/opening by reconstruction on `img`.

# Arguments
- `img::AbstractArray`: The input image.
- `se::AbstractArray`: The structuring element.
- `type::String`: The type of morphological operation to perform. Must be either "dilation" (close by reconstruction) or `"erosion"` (open by reconstruction).
- `invert::Bool=true`: Invert marker and mask before reconstruction.
"""
function reconstruct(img, se, type, invert::Bool=true)
    !(type == "dilation" || type == "erosion") &&
        throw(ArgumentError("Invalid type: $type. Must be 'dilation' or 'erosion'."))

    type == "dilation" && (
        morphed = to_uint8(
            IceFloeTracker.sk_morphology.dilation(img; footprint=collect(se))
        )
    )
    type == "erosion" && (
        morphed = to_uint8(IceFloeTracker.sk_morphology.erosion(img; footprint=collect(se)))
    )

    invert && (morphed = imcomplement(to_uint8(morphed)); img = imcomplement(img))

    type == "dilation" &&
        return IceFloeTracker.mreconstruct(IceFloeTracker.dilate, morphed, img)
    # dmw: the default for SK_Morphology is dilation, so even with the "erode" option, 
    # it's still doing reconstruction by dilation!
    # In SK morphology the default is a 3 by 3 box structuring element.
    return IceFloeTracker.sk_morphology.reconstruction(morphed, img)
end
