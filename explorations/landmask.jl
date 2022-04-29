using Images
using IceFloeTracker
using TiffImages
using LocalFilters
using ImageView
using ImageIO

tmp="hello"

test_image = load("/Users/tdivoll/ice-floe-tracker/existing_code/input/images/NE_Greenland.2020162.terra.250m.tif")
test_image2 = load("/Users/tdivoll/ice-floe-tracker/existing_code/input/images/NE_Greenland.2020163.terra.250m.tif")

foo = load("/Users/tdivoll/foo2.tiff")


land_mask = "/Users/tdivoll/ice-floe-tracker/existing_code/input/info/Land.tif"

mask_file = TiffImages.load(land_mask, mmap=true)

mask_file = dropdims(mask_file, dims = 3)

land_mask_binary = Gray.(mask_file) .== 0
land_mask_binary = LocalFilters.dilate(.!land_mask_binary, 50)
land_mask_binary = LocalFilters.closing(land_mask_binary, 15)

typeof(test_image)

test_dilation = (.!land_mask_binary .* test_image2)
save(test_dilation, "./dilated_test")

imshow(test_image)
using Colors

img = rand(2,2)

img2 = Gray.(img)

image(img)

channelview(test_dilation)

test = RGB.(test_dilation)

using TestImages
using ImageShow
image = testimage("mandrill")
imshow(image)

landmask_binary = IceFloeTracker.create_landmask("/Users/tdivoll/ice-floe-tracker/existing_code/input/info/Land.tif", 50, 15)
masked_images = IceFloeTracker.apply_landmask("/Users/tdivoll/ice-floe-tracker/existing_code/input/images/", landmask_binary)