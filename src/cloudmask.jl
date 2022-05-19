using Images
using ImageCore

test_image = load("test/data/NE_Greenland.2020162.aqua.250m.tiff")
sz = size(test_image)
top_test_image = test_image[ 1:2707, 1:4458]

cv = channelview(top_test_image)
# cv_int = trunc.(Int,cv.*255)
t = StackedView(zeroarray, cv[2,:,:], cv[3,:,:]) # ref_im3
colorview(RGB, t)

mask_aa = (cv[1,:,:] .< 0.784) # 200 on Int scale

mask_bb = (cv[2,:,:] .> 0.745) # 190 on Int scale

mask_cc = (mask_aa .&& mask_bb)

cv[1,:,:] = .!mask_cc .* cv[1,:,:]
cv[2,:,:] = .!mask_cc .* cv[2,:,:]

r = StackedView(cv[1,:,:], zeroarray, zeroarray)

mask_cloud_ice = Float64.(cv[1,:,:])./Float64.(cv[2,:,:])

mask_cloud_ice = (mask_cloud_ice .>= 0 .&& mask_cloud_ice .< 0.75)

ref_imclouds = (cv[1,:,:] .> 0.431)

ref_imclouds = .!mask_cloud_ice .* ref_imclouds

cv[1,:,:] = .!ref_imclouds .* cv[1,:,:]
cv[2,:,:] = .!ref_imclouds .* cv[2,:,:]
cv[3,:,:] = .!ref_imclouds .* cv[3,:,:]

ref_im7 = StackedView(zeroarray, cv[2,:,:], cv[3,:,:])

save("cloudmask_test.png", colorview(RGB, ref_im7))