"""
    create_cloudmask(reflectance_image)

Convert a 3-channel false color reflectance image to a 1-channel binary matrix; clouds = 0, else = 1.

# Arguments
- `reflectance_image`: corrected reflectance false color image - bands [7,2,1]

"""
function create_cloudmask(reflectance_image::Array{RGB{N0f8}})
  println("Setting thresholds")
  working_view = channelview(reflectance_image)
  orig_view = channelview(reflectance_image)
  orig_view_clouds = orig_view[1,:,:] .> 0.431 # intensity value 110
  mask_r = working_view[1,:,:] .< 0.784 # intensity value 200
  mask_g = working_view[2,:,:] .> 0.745 # intensity value 190
  println("Masking clouds and discriminatinhg cloud-ice")
  mask_rg = mask_r .&& mask_g
  working_view[1,:,:] = mask_rg .* working_view[1,:,:]
  working_view[2,:,:] = mask_rg .* working_view[2,:,:]
  cloud_ice = Float64.(working_view[1,:,:])./Float64.(working_view[2,:,:])
  mask_cloud_ice = cloud_ice .>= 0 .&& cloud_ice .< 0.75
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
function apply_cloudmask(reflectance_image::Array{RGB{N0f8}}, cloudmask::BitArray)
    masked_image = .!cloudmask .* reflectance_image
    image_view = channelview(masked_image)
    cloudmasked_image = StackedView(zeroarray, image_view[2,:,:], image_view[3,:,:])
    return cloudmasked_image
end