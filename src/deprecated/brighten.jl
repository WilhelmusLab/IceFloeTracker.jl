"""
    get_brighten_mask(equalized_gray_reconstructed_img, gamma_green)
# Arguments
- `equalized_gray_reconstructed_img`: The equalized gray reconstructed image (uint8 in Matlab).
- `gamma_green`: The gamma value for the green channel (also uint8).
# Returns
Difference equalized_gray_reconstructed_img - gamma_green clamped between 0 and 255.

# TODO: deprecate and do in place instead
"""
function get_brighten_mask(equalized_gray_reconstructed_img, gamma_green)
    return to_uint8(equalized_gray_reconstructed_img - gamma_green)
end

"""
    imbrighten(img, brighten_mask, bright_factor)
Brighten the image using a mask and a brightening factor.
# Arguments
- `img`: The input image.
- `brighten_mask`: A mask indicating the pixels to brighten.
- `bright_factor`: The factor by which to brighten the pixels.
# Returns
- The brightened image.

# TODO: deprecate and do in place instead
"""
function imbrighten(img, brighten_mask, bright_factor)
    img = Float64.(img)
    brighten_mask = brighten_mask .> 0
    img[brighten_mask] .= img[brighten_mask] * bright_factor
    return img = to_uint8(img)
end
