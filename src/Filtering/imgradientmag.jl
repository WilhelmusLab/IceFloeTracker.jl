import Images: imgradients
import Images.KernelFactors: sobel
"""
    imgradientmag(img)

Compute the gradient magnitude of an image using the Sobel operator.
"""
function imgradientmag(img, kernel=sobel)
    Gy, Gx = imgradients(img, kernel, "replicate")
    return hypot.(Gx, Gy)
end
