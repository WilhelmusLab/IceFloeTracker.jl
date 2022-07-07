"""
    normalize_image(truecolor_image, landmask, struct_elem; lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Adjusts truecolor land-masked image to highlight ice floe features. This function performs diffusion, adaptive histogram equalization, and sharpening, and returns a greyscale normalized image.

# Arguments
- `truecolor_image`: input image in truecolor
- `landmask`: bitmatrix landmask for region of interest
- `struct_elem`: structuring element for dilation
- `lambda`: speed of diffusion (0–0.25)
- `kappa`: conduction coefficient for diffusion (25–100)
- `niters`: number of iterations of diffusion
- `nbins`: number of bins during histogram equalization
- `rblocks`: number of row blocks to divide input image during equalization
- `cblocks`: number of column blocks to divide input image during equalization
- `clip`: threshold for clipping histogram bins (0–1); values closer to one minimize contrast enhancement, values closer to zero maximize contrast enhancement 
- `smoothing_param`: pixel radius for gaussian blurring (1–10)
- `intensity`: amount of sharpening to perform

"""
function normalize_image(
    truecolor_image::Matrix,
    landmask::BitMatrix,
    struct_elem::Matrix{Bool};
    lambda::Real=0.25,
    kappa::Real=90,
    niters::Int64=3,
    nbins::Int64=255,
    rblocks::Int64=8,
    cblocks::Int64=8,
    clip::Float64=0.95,
    smoothing_param::Int64=10,
    intensity::Float64=2.0,
)::Matrix
    gray_image = Float64.(Gray.(truecolor_image))
    image_diffused = diffusion(gray_image, 0.25, 75, 3)
    image_diffused_RGB = RGB.(image_diffused)
    masked_view = Float64.(channelview(image_diffused_RGB))

    image_equalized_1 = adjust_histogram(
        masked_view[1, :, :],
        AdaptiveEqualization(;
            nbins=255,
            rblocks=8,
            cblocks=8,
            minval=minimum(masked_view[1, :, :]),
            maxval=maximum(masked_view[1, :, :]),
            clip=0.8,
        ),
    )
    image_equalized_2 = adjust_histogram(
        masked_view[2, :, :],
        AdaptiveEqualization(;
            nbins=255,
            rblocks=8,
            cblocks=8,
            minval=minimum(masked_view[2, :, :]),
            maxval=maximum(masked_view[2, :, :]),
            clip=0.8,
        ),
    )
    image_equalized_3 = adjust_histogram(
        masked_view[3, :, :],
        AdaptiveEqualization(;
            nbins=255,
            rblocks=8,
            cblocks=8,
            minval=minimum(masked_view[3, :, :]),
            maxval=maximum(masked_view[3, :, :]),
            clip=0.8,
        ),
    )
    image_equalized = colorview(
        RGB, image_equalized_1, image_equalized_2, image_equalized_3
    )
    image_equalized_gray = Gray.(image_equalized)
    image_equalized_array = channelview(image_equalized_gray)

    image_smoothed = imfilter(image_equalized_gray, Kernel.gaussian(smoothing_param))
    image_smoothed_view = channelview(image_smoothed)

    image_sharpened =
        image_equalized_view .* (1 + intensity) .+ image_smoothed_view .* (-intensity)
    image_sharpened = max.(image_sharpened, 0.0)
    image_sharpened = min.(image_sharpened, 1.0)
    image_sharpened = colorview(Gray, image_sharpened)

    image_dilated = Images.dilate(image_sharpened, struct_elem)
    image_opened = Images.opening(complement.(image_dilated), complement.(image_sharpened))
    image_normalized_masked = IceFloeTracker.apply_landmask(image_opened, landmask)

    return image_normalized_masked
end
