# dmw: look for ways to avoid using these functions
function to_uint8(arr::AbstractMatrix{T}) where {T<:AbstractFloat}
    img = Int.(round.(arr, RoundNearestTiesAway))
    img = clamp.(img, 0, 255)
    return img
end

function to_uint8(arr::AbstractMatrix{T}) where {T<:Integer}
    img = clamp.(arr, 0, 255)
    return img
end

function to_uint8(num::T) where {T<:Union{AbstractFloat,Int,Signed}}
    num = Int(round(num, RoundNearestTiesAway))
    return clamp(num, 0, 255)
end

# dmw: use multiple dispatch, so that if the 2d function is called 
# with a 3D image we do this automatically
function anisotropic_diffusion_3D(I)
    rgbchannels = get_rgb_channels(I)

    for i in 1:3
        rgbchannels[:, :, i] .= anisotropic_diffusion_2D(rgbchannels[:, :, i])
    end

    return rgbchannels
end

# dmw: Move this to the image diffusion file
function anisotropic_diffusion_2D(
    # Implementation of the matlab 2D anisotropic diffusion filter default mode
    # by Carlos Paniagua
    I::AbstractMatrix{T}; gradient_threshold::Union{T,Nothing}=nothing, niter::Int=1
) where {T}
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
        diff_img_north =
            padded_img[0:(end - 1), 1:(end - 1)] .- padded_img[1:end, 1:(end - 1)]
        diff_img_east =
            padded_img[1:(end - 1), 1:end] .- padded_img[1:(end - 1), 0:(end - 1)]
        diff_img_nw = padded_img[0:(end - 2), 0:(end - 2)] .- I
        diff_img_ne = padded_img[0:(end - 2), 2:end] .- I
        diff_img_sw = padded_img[2:end, 0:(end - 2)] .- I
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
        flux_north_diff = flux_north[1:(end - 1), :] .- flux_north[2:end, :]
        flux_east_diff = flux_east[:, 2:end] .- flux_east[:, 1:(end - 1)]

        # Discrete PDE solution
        sum_ = (1 / (dd^2)) .* (flux_nw .+ flux_ne .+ flux_sw .+ flux_se)
        I = I .+ diffusion_rate .* (flux_north_diff .- flux_north_diff .+ sum_)
    end

    return I
end

# dmw: This function doesn't belong here
function imshow(img)
    if typeof(img) <: BitMatrix
        return Gray.(img)
    end
    return Gray.(img ./ 255)
end

function adapthisteq(img::Matrix{T}, nbins=256, clip=0.01) where {T}
    # Step 1: Normalize the image to [0, 1] based on its own min and max
    image_min, image_max = minimum(img), maximum(img)
    normalized_image = (img .- image_min) / (image_max - image_min)

    # Step 2: Apply adaptive histogram equalization. equalize_adapthist handles the tiling to 1/8 of the image size (equivalent to 8x8 blocks in MATLAB)
    equalized_image = sk_exposure.equalize_adapthist(
        normalized_image;
        clip_limit=clip,  # Equivalent to MATLAB's 'ClipLimit'
        nbins=nbins,         # Number of histogram bins. 255 is used to match the default in MATLAB script
    )

    # Step 3: Rescale the image back to the original range [image_min, image_max]
    final_image = sk_exposure.rescale_intensity(
        equalized_image; in_range="image", out_range=(image_min, image_max)
    )

    # Convert back to the original data type if necessary
    final_image = to_uint8(final_image)

    return final_image
end

"""
    get_rgb_channels(img)

Get the RBC (Red, Blue, and Green) channels of an image.

# Arguments
- `img`: The input image.

# Returns
An m x n x 3 array the Red, Blue, and Green channels of the input image.

"""
function get_rgb_channels(img)
    # TODO: might be able to use channelview instead
    # dmw: find ways to avoid casting to Int
    redc = red.(img) * 255
    greenc = green.(img) * 255
    bluec = blue.(img) * 255

    return cat(redc, greenc, bluec; dims=3)
end

"""
    rgb2gray(rgbchannels::Array{Float64, 3})

Convert an array of RGB channel data to grayscale in the range [0, 255].

Identical to MATLAB `rgb2gray` (https://www.mathworks.com/help/matlab/ref/rgb2gray.html).
"""
function rgb2gray(rgbchannels::Array{Float64,3}) # dmw: Can we set this up to return a Gray image instead of an int matrix?
    r, g, b = [to_uint8(rgbchannels[:, :, i]) for i in 1:3]
    # Reusing the r array to store the equalized gray image
    r .= to_uint8(0.2989 * r .+ 0.5870 * g .+ 0.1140 * b)
    return r
end

# dmw: let's figure out the difference between the three versions of this function
function _process_image_tiles(
    true_color_image,
    clouds_red,
    tiles,
    white_threshold,
    entropy_threshold,
    white_fraction_threshold,
)

    # Apply diffuse (anisotropic diffusion) to each channel of true color image
    true_color_diffused = IceFloeTracker.diffusion(float64.(true_color_image), 0.1, 75, 3)

    rgbchannels = get_rgb_channels(true_color_diffused)

    # For each tile, compute the entropy in the false color tile, and the fraction of white and black pixels
    for tile in tiles
        clouds_tile = clouds_red[tile...]
        entropy = Images.entropy(clouds_tile)
        whitefraction = sum(clouds_tile .> white_threshold) / length(clouds_tile)

        # If the entropy is above a threshold, and the fraction of white pixels is above a threshold, then apply histogram equalization to the tiles of each channel of the true color image. Otherwise, keep the original tiles.
        if entropy > entropy_threshold && whitefraction > white_fraction_threshold
            for i in 1:3
                eqhist = adapthisteq(rgbchannels[:, :, i][tile...])
                @view(rgbchannels[:, :, i])[tile...] .= eqhist
            end
        end
    end

    return rgbchannels
end

"""
    conditional_histeq(
    true_color_image,
    clouds_red,
    rblocks::Int,
    cblocks::Int,
    entropy_threshold::AbstractFloat=4.0,
    white_threshold::AbstractFloat=25.5,
    white_fraction_threshold::AbstractFloat=0.4,
)

Performs conditional histogram equalization on a true color image.

# Arguments
- `true_color_image`: The true color image to be equalized.
- `clouds_red`: The land/cloud masked red channel of the false color image.
- `rblocks`: The number of row-blocks to divide the image into for histogram equalization. Default is 8.
- `cblocks`: The number of column-blocks to divide the image into for histogram equalization. Default is 6.
- `entropy_threshold`: The entropy threshold used to determine if a block should be equalized. Default is 4.0.
- `white_threshold`: The white threshold used to determine if a pixel should be considered white. Default is 25.5.
- `white_fraction_threshold`: The white fraction threshold used to determine if a block should be equalized. Default is 0.4.

# Returns
The equalized true color image.

"""
function conditional_histeq(
    true_color_image,
    clouds_red,
    rblocks::Int,
    cblocks::Int,
    entropy_threshold::AbstractFloat=4.0,
    white_threshold::AbstractFloat=25.5,
    white_fraction_threshold::AbstractFloat=0.4,
)
    tiles = get_tiles(true_color_image; rblocks=rblocks, cblocks=cblocks)
    rgbchannels_equalized = _process_image_tiles(
        true_color_image,
        clouds_red,
        tiles,
        white_threshold,
        entropy_threshold,
        white_fraction_threshold,
    )

    return rgbchannels_equalized
end

"""
    conditional_histeq(
    true_color_image,
    clouds_red,
    side_length::Int,
    entropy_threshold::AbstractFloat=4.0,
    white_threshold::AbstractFloat=25.5,
    white_fraction_threshold::AbstractFloat=0.4,

)

Performs conditional histogram equalization on a true color image using tiles of approximately sidelength size `side_length`. If a perfect tiling is not possible, the tiling on the egde of the image is adjusted to ensure that the tiles are as close to `side_length` as possible. See `get_tiles(array, side_length)` for more details.
"""
function conditional_histeq(
    true_color_image,
    clouds_red,
    side_length::Int,
    entropy_threshold::AbstractFloat=4.0,
    white_threshold::AbstractFloat=25.5,
    white_fraction_threshold::AbstractFloat=0.4,
)
    side_length = IceFloeTracker.get_optimal_tile_size(side_length, size(true_color_image))

    tiles = IceFloeTracker.get_tiles(true_color_image, side_length)

    rgbchannels_equalized = _process_image_tiles(
        true_color_image,
        clouds_red,
        tiles,
        white_threshold,
        entropy_threshold,
        white_fraction_threshold,
    )

    return rgbchannels_equalized
end

function _get_false_color_cloudmasked(;
    false_color_image,
    prelim_threshold=110.0,
    band_7_threshold=200.0,
    band_2_threshold=190.0,
)
    mask_cloud_ice, clouds_view = IceFloeTracker._get_masks(
        false_color_image;
        prelim_threshold=prelim_threshold/255.,
        band_7_threshold=band_7_threshold/255.,
        band_2_threshold=band_2_threshold/255.,
        ratio_lower=0.0,
        ratio_upper=0.75
    )

    clouds_view[mask_cloud_ice] .= 0

    # remove clouds and land from each channel
    channels = Int.(channelview(false_color_image) * 255)

    # Apply the mask to each channel
    for i in 1:3
        @views channels[i, :, :][clouds_view] .= 0
    end

    return channels
end

"""
    rgb2gray(img::Matrix{RGB{Float64}})

Convert an RGB image to grayscale in the range [0, 255].
"""
function rgb2gray(img::Matrix{RGB{Float64}})
    return round.(Int, Gray.(img) * 255)
end

"""
    histeq(img)
    histeq(img; nbins=64)

Histogram equalization of `img` using `nbins` bins.
"""
function histeq(img::S; nbins=64)::S where {S<:AbstractArray{<:Integer}}
    return to_uint8(sk_exposure.equalize_hist(img; nbins=nbins) * 255)
end

function _imhist(img, rng)
    d = Dict(k => 0 for k in rng)
    for i in img
        d[i] = d[i] + 1
    end
    k, heights = collect.([keys(d), values(d)])
    order = sortperm(k)
    k, heights = k[order], heights[order]
    return k, heights
end

"""
    imhist(img, imgtype::AbstractString="uint8")

Compute the histogram of an image where each possible value is represented in the histogram. The function returns a tuple with the bins and counts of each bin.

# Example
```jldoctest; setup = :(using IceFloeTracker)
julia> img = [
    4 4 4 4 4
    3 4 5 4 3
    3 5 5 5 3
    3 4 5 4 3
    4 4 4 4 4
]

julia> bins, heights = imhist(img);

julia> [bins[heights .> 0] heights[heights .>0]] # display only non-zero bins and heights
3×2 Matrix{Int64}:
 3   6
 4  14
 5   5

"""
function imhist(img, imgtype::AbstractString="uint8")

    # TODO: add validation for arr: either uint8 0:255 or grayscale 0:1
    rng = imgtype == "uint8" ? range(0, 255) : range(0; stop=1, length=256)

    return _imhist(img, rng)
end
