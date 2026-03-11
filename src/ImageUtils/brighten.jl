import Images: AbstractGray, AbstractRGB, TransparentGray, TransparentColor
"""
    get_brighten_mask(equalized_gray_reconstructed_img, gamma_green)
# Arguments
- `equalized_gray_reconstructed_img`: The equalized gray reconstructed image (uint8 in Matlab).
- `gamma_green`: The gamma value for the green channel (also uint8).
# Returns
Difference equalized_gray_reconstructed_img - gamma_green clamped between 0 and 255.
"""
function get_brighten_mask(equalized_gray_reconstructed_img, gamma_green)
    return to_uint8(equalized_gray_reconstructed_img - gamma_green)
end

"""
    imbrighten(img, brighten_mask, bright_factor)
    imbrighten(img::AbstractArray{<:Union{AbstractRGB, TransparentColor, AbstractGray}},
     brighten_mask::Matrix{Bool}, bright_factor::Number))

Adjust image intensity within a masked region by multiplication with `bright_factor`. Despite the 
name, the function can also be used to selectively darken regions by supplying a bright factor between
0 and 1.

## Arguments
- `img`: The input image.
- `brighten_mask`: A mask indicating the pixels to brighten.
- `bright_factor`: The factor by which to brighten the pixels.
## Returns
- The brightened image.
"""
function imbrighten(img::AbstractArray{Int64}, brighten_mask, bright_factor)
    img = Float64.(img)
    brighten_mask = brighten_mask .> 0
    img[brighten_mask] .= img[brighten_mask] * bright_factor
    return img = to_uint8(img)
end

function imbrighten(
    img::AbstractArray{<:Union{AbstractRGB, TransparentColor, AbstractGray}},
    brighten_mask::AbstractArray{Bool}, bright_factor::Number
)
    _img = float64.(img)
    _img[brighten_mask] .= _img[brighten_mask] * bright_factor
    clamp01nan!(_img)
    return convert.(eltype(img), _img)
end
