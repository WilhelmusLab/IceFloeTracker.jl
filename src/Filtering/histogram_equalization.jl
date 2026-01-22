import Images: Images, RGB, float64, Gray, red, green, blue, adjust_histogram, AdaptiveEqualization
import ..skimage: sk_exposure
import ..ImageUtils: to_uint8, apply_to_channels

# dmw: use multiple dispatch, so that if the 2d function is called it automatically does this
function adapthisteq(img::Matrix{T}, nbins=256, clip=0.01) where {T}
    # TODO: Check whether the result is any different than using LinearStretching.
    # Step 1: Normalize the image to [0, 1] based on its own min and max
    image_min, image_max = minimum(img), maximum(img)
    normalized_image = (img .- image_min) / (image_max - image_min)

    # TODO: Consider setting up an AdaptiveEqualizationSKI that wraps the scikit image
    # adaptive equalization method and lets it be an argument to adjust_histogram.

    # Step 2: Apply adaptive histogram equalization. equalize_adapthist handles the tiling to 1/8 of the image size (equivalent to 8x8 blocks in MATLAB)
    equalized_image = sk_exposure.equalize_adapthist(
        normalized_image;
        clip_limit=clip,  # Equivalent to MATLAB's 'ClipLimit'
        nbins=nbins,         # Number of histogram bins. 255 is used to match the default in MATLAB script
    )

    # TODO: Check whether the result is any different than using LinearStretching.
    # Step 3: Rescale the image back to the original range [image_min, image_max]
    final_image = adjust_histogram(
        equalized_image, LinearStretching(nothing => (image_min, image_max))
    )

    # TODO: Let image types stay normed rather than going to UINT8.
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
    # TODO: Use channelview instead and don't cast to int.
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
    # TODO: Deprecate function and use Gray.() instead.
    r, g, b = [to_uint8(rgbchannels[:, :, i]) for i in 1:3]
    # Reusing the r array to store the equalized gray image
    r .= to_uint8(0.2989 * r .+ 0.5870 * g .+ 0.1140 * b)
    return r
end

# dmw: let's figure out the difference between the three versions of this function
"""
    conditional_histeq(
        image,
        clouds_red,
        tiles;
        entropy_threshold::Real=4.0,
        white_threshold::Real=25.5,
        white_fraction_threshold::Real=0.4,
    )

Performs conditional histogram equalization on a true color image.

# Arguments
- `image`: The true color image to be equalized.
- `clouds_red`: The land/cloud masked red channel of the false color image.
- `tiles`: the output from `get_tiles(image)` specifying the tiling to use on the image.
- `entropy_threshold`: The entropy threshold used to determine if a block should be equalized. Default is 4.0.
- `white_threshold`: The white threshold used to determine if a pixel should be considered white. Default is 25.5.
- `white_fraction_threshold`: The white fraction threshold used to determine if a block should be equalized. Default is 0.4.

# Returns
The equalized true color image.

"""
# TODO: Make this a "preprocessing workflow" function or functor.
function conditional_histeq(
    true_color_image,
    clouds_red,
    tiles;
    entropy_threshold::Real=4.0,
    white_threshold::Real=25.5,
    white_fraction_threshold::Real=0.4,
)

    # Apply Perona-Malik diffusion to each channel of true color image 
    # using the default inverse quadratic flux coefficient function
    pmd = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
    true_color_diffused = nonlinear_diffusion(float64.(true_color_image), pmd)

    # TODO: Replace the get_rgb_channels method with apply_to_channels
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

"""
    _channelwise_adapthisteq

    Broadcast adaptive histogram equalization to color channels of an image. Uses the 
    minimum and maximum value of each channel band instead of the 0 to 1 range from
    the default AdaptiveEqualization algorithm.

"""
function _channelwise_adapthisteq(img;
    nbins=256,
    rblocks=10,
    cblocks=10,
    clip=0.86)

    f_eq(img_channel) = adjust_histogram(img_channel,
        AdaptiveEqualization(;
            nbins=nbins,
            rblocks=rblocks,
            cblocks=cblocks,
            minval=minimum(img_channel),
            maxval=maximum(img_channel),
            clip=clip
            )
        )
    return apply_to_channels(img, f_eq) 
end