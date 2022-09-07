using Pkg;
Pkg.activate("explorations"); # delete later
using Revise

using IceFloeTracker
using Images
test_data_dir = "./test/test_inputs"
sharpened_test_file = "$(test_data_dir)/sharpened_test_image.png"
sharpened_image = load(sharpened_test_file)

count(x->((0.0001)<=x<=(0.4)), sharpened_image)

smallvals = .!(sharpened_image .< 0.4)

sharpened = (sharpened_image .* smallvals)

sharpened_int = trunc.(UInt8, sharpened.*255)

imgadj = adjust_gamma( sharpened, 10)

sharpened_cloudmasked = IceFloeTracker.apply_cloudmask(imgadj, cloudmask)


## skip .* 0.3, makes image too dark!


