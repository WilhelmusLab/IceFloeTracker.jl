"""
    normalize_image(truecolor_image; lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Adjusts truecolor land-masked image to highlight ice floe features. This function performs diffusion, adaptive histogram equalization, and sharpening, and returns a greyscale normalized image.

# Arguments
- `truecolor_image`: input image in truecolor
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
function normalize_image(truecolor_image::Matrix; lambda::Real=0.25, kappa::Real=90, niters::Int64=3, nbins::Int64=255, rblocks::Int64=8, cblocks::Int64=8, clip::Float64=0.95, smoothing_param::Int64=10, intensity::Float64=2.0)::Matrix
   
  test_data_dir = "../test/data"

  landmask = load("$(test_data_dir)/current_landmask.png")
  landmask_bm = convert(BitMatrix, landmask)

  gray_image = Float64.(Gray.(truecolor_image)) 

  imgdiffused = diffusion(gray_image, 0.25, 75, 3)

  imgdiffusedRGB = RGB.(imgdiffused)
  
  masked_v = Float64.(channelview(imgdiffusedRGB))
  
  imgeq_1 = adjust_histogram(masked_v[1,:,:], AdaptiveEqualization(nbins = 255, rblocks=8, cblocks=8, minval=minimum(masked_v[1,:,:]), maxval=maximum(masked_v[1,:,:]), clip=0.8))

  imgeq_2 = adjust_histogram(masked_v[2,:,:], AdaptiveEqualization(nbins = 255, rblocks=8, cblocks=8, minval=minimum(masked_v[2,:,:]), maxval=maximum(masked_v[2,:,:]), clip=0.8))

  imgeq_3 = adjust_histogram(masked_v[3,:,:], AdaptiveEqualization(nbins = 255, rblocks=8, cblocks=8, minval=minimum(masked_v[3,:,:]), maxval=maximum(masked_v[3,:,:]), clip=0.8))

  imgequalized = colorview(RGB, imgeq_1, imgeq_2, imgeq_3)

  img_equalized_gray = Gray.(imgequalized)
  img_smoothed = imfilter(img_equalized_gray, Kernel.gaussian(smoothing_param))
  img_equalized_array = channelview(img_equalized_gray)
  img_smoothed_array = channelview(img_smoothed)
  img_sharpened = img_equalized_array .* (1 + intensity) .+ img_smoothed_array .* (-intensity)
  img_sharpened = max.(img_sharpened, 0.0)
  img_sharpened = min.(img_sharpened, 1.0)
  img_sharpened = colorview(Gray, img_sharpened)
  
  strel_file2 = "$(test_data_dir)/se2.csv"
  struct_elem2 = readdlm(strel_file2, ',', Bool)
  
  img_dilated = Images.dilate(img_sharpened, struct_elem2)
  img_opened = Images.opening(complement.(img_dilated), complement.(img_sharpened))
  img_normalized_masked = IceFloeTracker.apply_landmask(img_opened, landmask_bm)
  save("test_output.png", img_normalized_masked)
  return img_normalized_masked
end

