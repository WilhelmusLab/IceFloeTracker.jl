"""
    create_cloudmask(reflectance_image)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix; clouds = 0, else = 1. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `ref_image`: corrected reflectance false color image - bands [7,2,1]

"""
function create_cloudmask(ref_image::Matrix{RGB{N0f8}}; prelim_threshold::Float64=0.431, red_threshold::Float64=0.784, green_threshold::Float64=0.745, ratio_lower::Float64=0.0, ratio_upper::Float64=0.75)
  println("Setting thresholds")
  ref_working_view = channelview(ref_image)
  ref_orig_view = channelview(ref_image)
  orig_view_clouds = ref_orig_view[1,:,:] .> prelim_threshold # intensity value 110
  mask_r = ref_working_view[1,:,:] .< red_threshold # intensity value 200
  mask_g = ref_working_view[2,:,:] .> green_threshold # intensity value 190
  # First find all the pixels that meet threshold logic in red and green channels
  println("Masking clouds and discriminatinhg cloud-ice")
  # Next find pixels that meet both thresholds and mask them from red and green channels
  mask_rg = mask_r .&& mask_g
  ref_working_view[1,:,:] = mask_rg .* ref_working_view[1,:,:]
  ref_working_view[2,:,:] = mask_rg .* ref_working_view[2,:,:]
  cloud_ice = Float64.(ref_working_view[1,:,:])./Float64.(ref_working_view[2,:,:])
  mask_cloud_ice = cloud_ice .>= ratio_lower .&& cloud_ice .< ratio_upper
  println("Creating final cloudmask")
  cloudmask = .!mask_cloud_ice .* orig_view_clouds
  return cloudmask
end

"""
    apply_cloudmask(reflectance_image, cloudmask)

Zero out pixels containing clouds where clouds and ice are not discernable.

# Arguments
- `reflectance_image`: corrected reflectance false color image - bands [7,2,1]
- `cloudmask`: binary cloudmask with clouds = 0, else = 1

"""
function apply_cloudmask(ref_image::Matrix{RGB{N0f8}}, cloudmask::BitMatrix)::Matrix{RGB{N0f8}}
    masked_image = .!cloudmask .* ref_image
    image_view = channelview(masked_image)
    cloudmasked_view = StackedView(zeroarray, image_view[2,:,:], image_view[3,:,:])
    cloudmasked_image = colorview(RGB, cloudmasked_view)
    return cloudmasked_image
end