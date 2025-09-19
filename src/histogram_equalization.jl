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
function rgb2gray(rgbchannels::Array{Float64,3}) 
# dmw: Can we set this up to return a Gray image instead of an int matrix?
# dmw: Check whether the coefficients differ substantially -- we can make it an alternative to the Gray function, like GrayMLB and apply it like Gray.()
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

    # Apply Perona-Malik diffusion to each channel of true color image 
    # using the default inverse quadratic flux coefficient function
    pmd = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
    true_color_diffused = IceFloeTracker.nonlinear_diffusion(float64.(true_color_image), pmd)

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

"""
Private function for testing the conditional adaptive histogram equalization workflow.
"""
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
        ratio_offset=0.0,
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
