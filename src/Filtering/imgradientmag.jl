# import Images: centered
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

# function imgradientmag(img)
#     h = centered([-1 0 1; -2 0 2; -1 0 1]')
#     Gx_future = Threads.@spawn imfilter(img, h', "replicate")
#     Gy_future = Threads.@spawn imfilter(img, h, "replicate")
#     Gx = fetch(Gx_future)
#     Gy = fetch(Gy_future)
#     return hypot.(Gx, Gy)
# end
