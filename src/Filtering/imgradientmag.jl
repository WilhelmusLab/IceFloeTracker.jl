import Images: imgradients
import Images.KernelFactors: sobel
"""
    imgradientmag(img, kernel=Images.KernelFactors.sobel)

Compute the gradient magnitude of an image using the specified operator. Wrapper for ImageFiltering `imgradients`.

"""
function imgradientmag(img, kernel=sobel)
    Gy, Gx = imgradients(img, kernel, "replicate")
    return hypot.(Gx, Gy)
end
