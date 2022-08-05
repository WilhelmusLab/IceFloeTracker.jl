using Pkg;
Pkg.activate("explorations"); # delete later
using Revise # delete later
using Images
using Clustering
using ImageView
using IceFloeTracker
using ImageSegmentation
using Peaks
using StatsBase

test_data_dir = "./test/test_inputs"
ice_water_discrim_test_file = "$(test_data_dir)/ice_water_discrim_image.png"

ice_water_matlab = "$(test_data_dir)/matlab_ice_water_discrim.png"
mat_he = load(ice_water_matlab)
#input = ice_water_discriminated_image

he = load(ice_water_discrim_test_file)

# # Fuzzy C-means - probably works but how do we visualize?
# k = 4
# fuzziness = 1.01
# he = (RGB.(he))

# r = fuzzy_cmeans(he, k, fuzziness)
# #centers = colorview(RGB, r.centers)
# pic1 = Gray.(r.centers[1] * reshape(r.weights[:, 1], axes(he)))
# pic2 = N0f8.(Gray.(r.centers[1] * reshape(r.weights[:, 2], axes(he)))*1.5)
# pic3 = Gray.(r.centers[1] * reshape(r.weights[:, 3], axes(he)))

# IceFloeTracker.@persist pic2 "fuzzy_c.png"

## Felzenszwalb region growing
# segments = felzenszwalb(he, 45000, 10)

# function get_random_color(seed)
#     Random.seed!(seed)
#     rand(RGB{N0f8})
# end
# imshow(map(i->get_random_color(i), labels_map(segments)))

# # Fast-scanning
# fast = fast_scanning(he_1, 0.1)

# # unseeded region growing
# usg = unseeded_region_growing(he_1, 0.1)

# K-means
# make landmask white rather than black
# @. he = ifelse(he==0, 1, he)
he = Array{Float32,2}(he)
Nx, Ny = size(he)
dat = reshape(he, 1, Nx * Ny)
R = kmeans(dat, 4; maxiter=50, display=:iter, init=:kmpp)
a = assignments(R)
segmented = Gray.(((reshape(a, Nx, Ny)) .- 1) ./ 3) ## pixel_labels in matlab

## Make ice masks
reflectance_test_image_file = "$(test_data_dir)/NE_Greenland.2020162.aqua.250m.tiff"
current_landmask_file = "$(test_data_dir)/current_landmask.png"
landmask_bitmatrix = convert(BitMatrix, load(current_landmask_file))
ref_image = load(reflectance_test_image_file)[test_region...]
cv = channelview(ref_image)
mask_ice_1 = cv[1, :, :] .< 5 / 255
mask_ice_2 = cv[2, :, :] .> 230 / 255
mask_ice_3 = cv[3, :, :] .> 240 / 255

ice_labels = IceFloeTracker.apply_landmask(
    Gray.(mask_ice_1 .|| mask_ice_2 .|| mask_ice_3), landmask_bitmatrix
)

if isempty(ice_labels)
    mask_ice_1 = cv[1, :, :] .< 10 / 255
    mask_ice_3 = cv[3, :, :] .> 190 / 255
    ice_labels = IceFloeTracker.apply_landmask(
        Gray.(mask_ice_1 .|| mask_ice_2 .|| mask_ice_3), landmask_bitmatrix
    )
    if isempty(ice_labels)
        ref22 = cv[2, :, :]
        ref33 = cv[3, :, :]
        ref22[ref22 .< 75 / 255] .= 0
        ref33[ref33 .< 75 / 255] .= 0
        edges22, counts22 = ImageContrastAdjustment.build_histogram(ref22)
        locs22, pks22 = Peaks.findmaxima(counts22)
        locs22 = sort(locs22; rev=true)
        peak1 = locs22[2]
        edges33, counts33 = ImageContrastAdjustment.build_histogram(ref33)
        locs33, pks33 = Peaks.findmaxima(counts33)
        locs33 = sort(locs33; rev=true)
        peak2 = locs33[2]
        mask_ice_2 = cv[2, :, :] .> peak1 / 255
        mask_ice_3 = cv[3, :, :] .> peak2 / 255
        ice_labels = IceFloeTracker.apply_landmask(
            Gray.(mask_ice_1 .|| mask_ice_2 .|| mask_ice_3), landmask_bitmatrix
        )
        segmented_masked = segmented .* ice_labels
        nlabel = mean(segmented_masked)
    end
end
