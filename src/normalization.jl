"""
    normalize_image(landmasked_image; lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

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
function normalize_image(landmasked_image::Matrix; lambda::Real=0.25, kappa::Real=50, niters::Int64=3, nbins::Int64=255, rblocks::Int64=8, cblocks::Int64=8, clip::Float64=0.75, smoothing_param::Int64=10, intensity::Float64=2.0)::Matrix
  # println("Applying Nonlinear diffusion filtering") 
  test_data_dir = "../test/data"
  test_region = (1:2707, 1:4458)
  landmask_file = "$(test_data_dir)/landmask.tiff"
  lm_image = load(landmask_file)[test_region...]
  strel_file = "$(test_data_dir)/se.csv"
  num_pixels_closing = 50
  struct_elem = readdlm(strel_file, ',', Bool)
  landmask = IceFloeTracker.create_landmask(lm_image, struct_elem; num_pixels_closing=num_pixels_closing)

  masked_v = Float64.(channelview(landmasked_image))
  
  imgeq_1 = adjust_histogram(masked_v[1,:,:], AdaptiveEqualization(nbins = 255, rblocks=8, cblocks=8, minval=minimum(masked_v[1,:,:]), maxval=maximum(masked_v[1,:,:]), clip=0.75))

  imgeq_2 = adjust_histogram(masked_v[2,:,:], AdaptiveEqualization(nbins = 255, rblocks=8, cblocks=8, minval=minimum(masked_v[2,:,:]), maxval=maximum(masked_v[2,:,:]), clip=0.75))

  imgeq_3 = adjust_histogram(masked_v[3,:,:], AdaptiveEqualization(nbins = 255, rblocks=8, cblocks=8, minval=minimum(masked_v[3,:,:]), maxval=maximum(masked_v[3,:,:]), clip=0.75))

  imgequalized = colorview(RGB, imgeq_1, imgeq_2, imgeq_3)

  # println("Applying CLAHE - Contrast Limited Adaptive Histogram Equalization")
 
  # println("Sharpening Image")
  img_equalized_gray = Gray.(imgequalized)
  img_smoothed = imfilter(img_equalized_gray, Kernel.gaussian(smoothing_param))
  img_equalized_array = channelview(img_equalized_gray)
  img_smoothed_array = channelview(img_smoothed)
  img_sharpened = @. img_equalized_array * (1 + intensity) + img_smoothed_array * (-intensity)
  img_sharpened = max.(img_sharpened, 0.0)
  img_sharpened = min.(img_sharpened, 0.1)
  img_sharpened = colorview(Gray, img_sharpened)
  
  strel_file2 = """$(test_data_dir)/se2.csv"""
  struct_elem2 = readdlm(strel_file2, ',', Bool)
  
  Iobrd = Images.dilate(img_sharpened, struct_elem2)
  Iobrcbr = Images.opening(complement.(Iobrd), complement.(img_sharpened))
  Iobrcbr_masked = IceFloeTracker.apply_landmask(Iobrcbr, landmask)
  return Iobrcbr_masked
end

