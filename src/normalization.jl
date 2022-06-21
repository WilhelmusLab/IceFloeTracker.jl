"""
    normalize_image(landmasked_image, lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Adjusts truecolor land-masked image to highlight ice floe features. This function performs diffusion, adaptive histogram equalization, and sharpening, and returns a greyscale normalized image.

# Arguments
- `landmasked_image`: land-masked image in truecolor
- `lambda`: speed of diffusion (0–0.25)
- `kappa`: conduction coefficient for diffusion (25–100)
- `niters`: number of iterations of diffusion
- `nbins`: number of bins during histogram equalization
- `rblocks`: number of row blocks to divide input image during equalization
- `cblocks`: number of column blocks to divide input image during equalization
- `clip`: threshold for clipping histogram bins (0–1); values closer to one minimize contrast enhancement, values closer to zero maximize contrast enhancement 
- `smoothing_param`: pixel radius for gaussian blurring (1–10)
- `intensity`: amount of sharpening to perform

"""
function normalize_image(landmasked_image::Matrix{RGB{N0f8}}; lambda::Real=0.25, kappa::Real=50, niters::Int64=3, nbins::Int64=255, rblocks::Int64=9, cblocks::Int64=6, clip::Float64=0.75, smoothing_param::Int64=10, intensity::Float64=2.0)::Matrix{Gray{N0f8}}

  println("Applying Nonlinear diffusion filtering") 
  masked_view = Float64.(channelview(landmasked_image))
  include("../src/diffusion.jl");
  ch1 = diffusion((masked_view[1,:,:]), lambda, kappa, niters);
  ch2 = diffusion((masked_view[2,:,:]), lambda, kappa, niters);
  ch3 = diffusion((masked_view[3,:,:]), lambda, kappa, niters);
  diffused_image = colorview(RGB, ch1, ch2, ch3)

  println("Applying CLAHE - Contrast Limited Adaptive Histogram Equalization")
 
  img_eq = Gray.(adjust_histogram(diffused_image, AdaptiveEqualization(nbins=nbins, rblocks=rblocks, cblocks=cblocks, clip=clip)))
 
  println("Sharpening Image")
  
  img_b = imfilter(img_eq, Kernel.gaussian(smoothing_param))
  img_eq_array = chanelview(img_eq)
  img_b_array = chanelview(img_b)
  img_sharp = @. img_eq_array * (1 + intensity) + img_b_array * (-intensity)
  img_sharp = max.(img_sharp, 0.0)
  img_sharp = min.(img_sharp, 0.1)
  img_sharp = colorview(Gray, img_sharp)
  return img_sharp
end

