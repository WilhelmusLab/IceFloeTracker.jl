using StatsBase
using TiledIteration
using DelimitedFiles
using IceFloeTracker
using Images

function shannon_entropy(img)
    # Compute the histogram and normalize
    hist = fit(Histogram, vec(img), 0:255)
    p = hist.weights / sum(hist.weights)

    # Filter out zero probabilities and compute entropy
    p = p[p.>0]
    H = -sum(p .* log2.(p))

    return H
end

function anisotropic_diffusion_3D(I)
    rgbchannels = getrbc_channels(I)

    for i in 1:3
        rgbchannels[:, :, i] .= anisotropic_diffusion_2D(rgbchannels[:, :, i])
    end

    return rgbchannels

end

function anisotropic_diffusion_2D(I::AbstractMatrix{T}; gradient_threshold::Union{T,Nothing}=nothing, niter::Int=1) where {T}
    if eltype(I) <: Int
        I = Gray.(I ./ 255)
    end

    # Determine the gradient threshold if not provided
    if gradient_threshold === nothing
        dynamic_range = maximum(I) - minimum(I)
        gradient_threshold = 0.1 * dynamic_range
    end

    # Padding the image (corrected)
    padded_img = padarray(I, Pad(:replicate, (1, 1)))
    dd = sqrt(2)
    diffusion_rate = 1 / 8  # Fixed for maximal connectivity (8 neighbors)

    for _ in 1:niter
        # These are zero-indexed offset arrays
        diff_img_north = padded_img[0:end-1, 1:end-1] .- padded_img[1:end, 1:end-1]
        diff_img_east = padded_img[1:end-1, 1:end] .- padded_img[1:end-1, 0:end-1]
        diff_img_nw = padded_img[0:end-2, 0:end-2] .- I
        diff_img_ne = padded_img[0:end-2, 2:end] .- I
        diff_img_sw = padded_img[2:end, 0:end-2] .- I
        diff_img_se = padded_img[2:end, 2:end] .- I

        # Exponential conduction coefficients
        conduct_coeff_north = exp.(-(abs.(diff_img_north) ./ gradient_threshold) .^ 2)
        conduct_coeff_east = exp.(-(abs.(diff_img_east) ./ gradient_threshold) .^ 2)
        conduct_coeff_nw = exp.(-(abs.(diff_img_nw) ./ gradient_threshold) .^ 2)
        conduct_coeff_ne = exp.(-(abs.(diff_img_ne) ./ gradient_threshold) .^ 2)
        conduct_coeff_sw = exp.(-(abs.(diff_img_sw) ./ gradient_threshold) .^ 2)
        conduct_coeff_se = exp.(-(abs.(diff_img_se) ./ gradient_threshold) .^ 2)


        # Flux calculations
        flux_north = conduct_coeff_north .* diff_img_north
        flux_east = conduct_coeff_east .* diff_img_east
        flux_nw = conduct_coeff_nw .* diff_img_nw
        flux_ne = conduct_coeff_ne .* diff_img_ne
        flux_sw = conduct_coeff_sw .* diff_img_sw
        flux_se = conduct_coeff_se .* diff_img_se

        # Back to regular 1-indexed arrays
        flux_north_diff = flux_north[1:end-1, :] .- flux_north[2:end, :]
        flux_east_diff = flux_east[:, 2:end] .- flux_east[:, 1:end-1]

        # Discrete PDE solution
        sum_ = (1 / (dd^2)) .* (flux_nw .+ flux_ne .+ flux_sw .+ flux_se)
        I = I .+ diffusion_rate .* (flux_north_diff .- flux_north_diff .+ sum_)

    end

    return I
end


function imshow(img)
    if typeof(img) <: BitMatrix
        return Gray.(img)
    end
    Gray.(img ./ 255)
end

function to_uint8(img)
    img = UInt8.(round.(img))
    img = clamp.(img, 0, 255)
    return img
end

function adapthisteq(img::Matrix{T}, nbins=256, clip=0.01) where {T}
    # Step 1: Normalize the image to [0, 1] based on its own min and max
    image_min, image_max = minimum(img), maximum(img)
    normalized_image = (img .- image_min) / (image_max - image_min)

    # Step 2: Apply adaptive histogram equalization. equalize_adapthist handles the tiling to 1/8 of the image size (equivalent to 8x8 blocks in MATLAB)
    equalized_image = sk_exposure.equalize_adapthist(
        normalized_image,
        clip_limit=clip,  # Equivalent to MATLAB's 'ClipLimit'
        nbins=nbins         # Number of histogram bins. 255 is used to match the default in MATLAB script
    )

    # Step 3: Rescale the image back to the original range [image_min, image_max]
    final_image = sk_exposure.rescale_intensity(equalized_image, in_range="image", out_range=(image_min, image_max))

    # Convert back to the original data type if necessary
    final_image = to_uint8(final_image)

    return final_image
end

"""
    getrbc_channels(img)

Get the RBC (Red, Blue, and Green) channels of an image.

# Arguments
- `img`: The input image.

# Returns
An m x n x 3 array the Red, Blue, and Green channels of the input image.

"""
function getrbc_channels(img)
    # TODO: might be able to use channelview instead
    redc = red.(img) * 255
    greenc = green.(img) * 255
    bluec = blue.(img) * 255

    return cat(redc, greenc, bluec, dims=3)
end


"""
    conditional_histeq(true_color_img, clouds, landmask, rblocks::Int=8, cblocks::Int=6, entropy_threshold::AbstractFloat=4.0, white_threshold::AbstractFloat=25.5, white_fraction_threshold::AbstractFloat=0.4)

Performs conditional histogram equalization on a true color image.

# Arguments
- `true_color_img`: The true color image to be equalized.
- `clouds_red`: Preproccesed red channel of false color (clouds) image for this purpose. Used to determine which regions of the true color image should be equalized.
- `landmask`: The land mask indicating the land regions in the image.
- `rblocks`: The number of row-blocks to divide the image into for histogram equalization. Default is 8.
- `cblocks`: The number of column-blocks to divide the image into for histogram equalization. Default is 6.
- `entropy_threshold`: The entropy threshold used to determine if a block should be equalized. Default is 4.0.
- `white_threshold`: The white threshold used to determine if a pixel should be considered white. Default is 25.5.
- `white_fraction_threshold`: The white fraction threshold used to determine if a block should be equalized. Default is 0.4.

# Returns
The equalized true color image.

"""
function conditional_histeq(
    true_color_img,
    clouds_red,
    landmask,
    rblocks::Int=8,
    cblocks::Int=6,
    entropy_threshold::AbstractFloat=4.0,
    white_threshold::AbstractFloat=25.5,
    white_fraction_threshold::AbstractFloat=0.4)

    # 1. Apply diffuse (anisotropic diffusion) to each channel of true color image
    # TODO: See about using 8-point connectivity as in the original MATLAB script
    true_color_diffused = IceFloeTracker.diffusion(float64.(true_color_img), 0.1, 75, 3)

    # 2. Apply landmask to clouds/falscolor image
    IceFloeTracker.apply_landmask!(clouds_red, landmask.dilated)

    # 3. Get tiles
    # TODO: Define tile size based on the desired number of pixels per tile (WIP)
    rtile, ctile = size(clouds_red)
    tile_size = Tuple{Int,Int}((rtile / rblocks, ctile / cblocks))
    tiles = TileIterator(axes(clouds_red), tile_size)

    rgbchannels = getrbc_channels(true_color_diffused)

    # 4. For each tile, compute the entropy in the falscolor tile, and the fraction of white and black pixels
    for tile in tiles
        clouds_tile = clouds_red[tile...]
        entropy = shannon_entropy(clouds_tile)
        whitefraction = sum(clouds_tile .> white_threshold) / length(clouds_tile)

        # 5. If the entropy is above a threshold, and the fraction of white pixels is above a threshold, then apply histogram equalization to the tiles of each channel of the true color image. Otherwise, keep the original tiles.
        if entropy > entropy_threshold && whitefraction > white_fraction_threshold
            for i in 1:3
                eqhist = adapthisteq(rgbchannels[:, :, i][tile...])
                @view(rgbchannels[:, :, i])[tile...] .= eqhist
            end
        end

    end

    return rgbchannels
end