using IceFloeTracker
using Images
img = rand(RGBA{N0f8}, 2, 2)
@show img
@show IceFloeTracker.diffusion(img, 0.1, 75, 3)