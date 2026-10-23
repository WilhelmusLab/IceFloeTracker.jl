import Images:
    RGB, float64, Gray, entropy, red, green, blue, adjust_histogram, AdaptiveEqualization
import ..ImageUtils: to_uint8, apply_to_channels, @tilewise

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

    rgbchannels = get_rgb_channels(true_color_diffused)

    return @tilewise tiles _conditional_histeq_tile(
        $rgbchannels,
        $clouds_red;
        entropy_threshold=entropy_threshold,
        white_threshold=white_threshold,
        white_fraction_threshold=white_fraction_threshold,
    ) into = rgbchannels
end

"""
    _conditional_histeq_tile(rgb, clouds; entropy_threshold, white_threshold, white_fraction_threshold)

Equalize each channel of the `h×w×3` tile `rgb` with contrast-limited adaptive histogram
equalization if the false color tile `clouds` has entropy above `entropy_threshold` and a fraction
of pixels above `white_threshold` greater than `white_fraction_threshold`. Otherwise return `rgb`.
"""
function _conditional_histeq_tile(
    rgb, clouds; entropy_threshold, white_threshold, white_fraction_threshold
)
    whitefraction = sum(clouds .> white_threshold) / length(clouds)
    (entropy(clouds) > entropy_threshold && whitefraction > white_fraction_threshold) ||
        return rgb
    clahe = ContrastLimitedAdaptiveHistogramEqualization(; clip=2.0)
    return stack(adjust_histogram(rgb[:, :, i] / 255.0, clahe) .* 255 for i in 1:3)
end

"""
    rgb2gray(img::Matrix{RGB{Float64}})

Convert an RGB image to grayscale in the range [0, 255].
"""
function rgb2gray(img::Matrix{RGB{Float64}})
    return round.(Int, Gray.(img) * 255)
end
