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
- `clip`: tuple (one per channel) with thresholds for clipping histogram bins (0–1); values closer to one minimize contrast enhancement, values closer to zero maximize contrast enhancement 
- `smoothing_param`: pixel radius for gaussian blurring (1–10)
- `intensity`: amount of sharpening to perform

"""
function normalize_image(
    truecolor_image::Matrix{RGB{Float64}},
    landmask::BitMatrix,
    struct_elem::Matrix{Bool};
    lambda::Real=0.25,
    kappa::Real=75,
    niters::Int64=3,
    nbins::Int64=255,
    rblocks::Int64=8,
    cblocks::Int64=8,
    clip::Tuple{Float64, Float64, Float64}= (0.95, 0.8, 0.8),
    smoothing_param::Int64=10,
    intensity::Float64=2.0,
)::Tuple{Matrix{Gray{Float64}},Matrix{Gray{Float64}}}
    image_sharpened = imsharpen(truecolor_image, lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

    # Apply landmask and turn to gray
    image_sharpened_gray = imsharpen_gray(image_sharpened, landmask)

    image_dilated = ImageMorphology.dilate(image_sharpened_gray, struct_elem)

    image_opened = ImageMorphology.opening(
        complement.(image_dilated), complement.(image_sharpened)
    )
    image_normalized_masked = IceFloeTracker.apply_landmask(image_opened, landmask)

    return image_sharpened_gray, image_normalized_masked
end

function _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)
    adjust_histogram(
        masked_view,
        AdaptiveEqualization(;
            nbins=nbins,
            rblocks=rblocks,
            cblocks=cblocks,
            minval=minimum(masked_view),
            maxval=maximum(masked_view),
            clip=clip,
        ),
    )
end

function imsharpen(truecolor_image, lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)
    gray_image = Float64.(Gray.(truecolor_image))
    image_diffused = diffusion(gray_image, lambda, kappa, niters)
    image_diffused_RGB = RGB.(image_diffused)
    masked_view = Float64.(channelview(image_diffused_RGB))
      
    eq = [_adjust_histogram(masked_view[i,:,:],nbins, rblocks, cblocks, clip[i]) for i=1:3]
    image_equalized = colorview(RGB, eq...)
    image_equalized_gray = Gray.(image_equalized)
    image_equalized_view = channelview(image_equalized_gray)

    image_smoothed = imfilter(image_equalized_gray, Kernel.gaussian(smoothing_param))
    image_smoothed_view = channelview(image_smoothed)

    image_sharpened =
        image_equalized_view .* (1 + intensity) .+ image_smoothed_view .* (-intensity)
    image_sharpened = max.(image_sharpened, 0.0)
    image_sharpened = min.(image_sharpened, 1.0)
end

function imsharpen_gray(imgsharpened, landmask)
    image_sharpened_landmasked = apply_landmask(imgsharpened, landmask)
    colorview(Gray, image_sharpened_landmasked)
end
