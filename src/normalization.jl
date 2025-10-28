"""
    normalize_image(image_sharpened, image_sharpened_gray, landmask, struct_elem;)

Adjusts sharpened land-masked image to highlight ice floe features.

Does reconstruction and landmasking to `image_sharpened`.

# Arguments
- `image_sharpened`: sharpened image (output of `imsharpen`)
- `image_sharpened_gray`: grayscale, landmasked sharpened image (output of `imsharpen_gray(image_sharpened)`)
- `landmask`: landmask for region of interest
- `struct_elem`: structuring element for dilation

"""
function normalize_image(
    image_sharpened::Matrix{Float64},
    image_sharpened_gray::T,
    landmask::BitMatrix,
    struct_elem;
)::Matrix{Gray{Float64}} where {T<:AbstractMatrix{Gray{Float64}}}
    image_dilated = dilate(image_sharpened_gray, struct_elem)

    image_reconstructed = mreconstruct(
        dilate, complement.(image_dilated), complement.(image_sharpened)
    )
    return IceFloeTracker.apply_landmask(image_reconstructed, landmask)
end

function normalize_image(
    image_sharpened::Matrix{Float64},
    image_sharpened_gray::Matrix{Gray{Float64}},
    landmask::BitMatrix,
)::Matrix{Gray{Float64}}
    return normalize_image(
        image_sharpened, image_sharpened_gray, landmask, strel_diamond((5, 5))
    )
end

# TODO: Remove function, replace with direct use of landmask and colorview.
"""
    imsharpen_gray(imgsharpened, landmask)

Apply landmask and return Gray type image in colorview for normalization.

"""
function imsharpen_gray(
    imgsharpened::Matrix{Float64}, landmask::AbstractArray{Bool}
)::Matrix{Gray{Float64}}
    image_sharpened_landmasked = apply_landmask(imgsharpened, landmask)
    return colorview(Gray, image_sharpened_landmasked)
end

# TODO: Remove once the workflow is all normed images
function adjustgamma(img, gamma=1.5, asuint8=true)
    if maximum(img) > 1
        img = img ./ 255
    end

    adjusted = adjust_histogram(img, GammaCorrection(gamma))

    if asuint8
        adjusted = Int.(round.(adjusted * 255, RoundNearestTiesAway))
    end

    return adjusted
end

# TODO: Remove function
function imbinarize(img)
    f = AdaptiveThreshold(img) # infer the best `window_size` using `img`
    return binarize(img, f)
end
