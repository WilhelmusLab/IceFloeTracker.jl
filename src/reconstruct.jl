# WIP
function imcomplement(img::Matrix{T}) where T<:Union{Unsigned, Int}
    return 255 .- img
end

function imcomplement(img::Matrix{Gray{Float64}})
    return 1 .- img
end




function reconstruct(img, se, type, invert)

    if type == "dilation"
        morphed = IceFloeTracker.MorphSE.erode(img, se)
    elseif type == "erosion"
        morphed = IceFloeTracker.MorphSE.dilate(img, se)
    end

    if invert == 1
        morphed = imcomplement(morphed);
        img = imcomplement(img);
    end

    reconstructed = imreconstruct(morphed, img);

end