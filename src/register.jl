greaterthan05(x) = x .> 0.5 # used for the image resize step and for binarizing images
function imrotate_bin(x, r)
    return greaterthan05(collect(imrotate(x, r, axes(x); method=BSpline(Constant()))))
end
function imrotate_bin_nocrop(x, r)
    return greaterthan05(collect(imrotate(x, r; method=BSpline(Constant()))))
end
imrotate_bin_clockwise_radians(x, r) = imrotate_bin(x, r)
imrotate_bin_counterclockwise_radians(x, r) = imrotate_bin(x, -r)
imrotate_bin_clockwise_degrees(x, r) = imrotate_bin_clockwise_radians(x, deg2rad(r))
imrotate_bin_counterclockwise_degrees(x, r) = imrotate_bin_counterclockwise_radians(x, deg2rad(r))



# Functions used for the SD minimization
"""
Pad images by zeros based on the size of the larger of the two images.
"""
function pad_images(im1, im2)
    n = max(size(im1)..., size(im2)...)
    im_padded = [collect(padarray(im, Fill(0, (n, n), (n, n)))) for im in [im1, im2]]
    return im_padded[1], im_padded[2]
end

"""
Calculate the centroid of a binary image. If 'rounded', return the
nearest integer.
"""
function compute_centroid(im::AbstractArray{Bool}; rounded=false)
    xi = 0
    yi = 0
    R = sum(im .> 0)
    for idx in CartesianIndices(im)
        if im[idx] > 0
            ii, jj = Tuple(idx)
            xi += ii
            yi += jj
        end
    end

    x0, y0 = sum(xi) / R, sum(yi) / R
    if rounded
        return round(Int32, x0), round(Int32, y0)
    else
        return x0, y0
    end
end

"""
Align images by selecting and cropping so that r1, c1 and r2, c2 are the center.
These values are expected to be the (integer) centroid of the image. These images
should already be padded so that there is no danger of cutting into the floe shape.
"""
function crop_to_shared_centroid(im1, im2)
    r1, c1 = compute_centroid(im1; rounded=true)
    r2, c2 = compute_centroid(im2; rounded=true)

    n1, m1 = size(im1)
    n2, m2 = size(im2)
    new_halfn = minimum([minimum([r1, n1 - r1]), minimum([r2, n2 - r2])])
    new_halfm = minimum([minimum([c1, m1 - c1]), minimum([c2, m2 - c2])])

    # check notation: how does julia interpret start and end of array index?
    im1_cropped = im1[
        (1+r1-new_halfn):(r1+new_halfn), (1+c1-new_halfm):(c1+new_halfm)
    ]
    im2_cropped = im2[
        (1+r2-new_halfn):(r2+new_halfn), (1+c2-new_halfm):(c2+new_halfm)
    ]

    return im1_cropped, im2_cropped
end

"""
Computes the shape difference between im_reference and im_target for each angle in test_angles.
The reference image is held constant, while the target image is rotated. The test_angles are interpreted
as the angle of rotation from target to reference, so to find the best match, we rotate the reverse
direction. A perfect match at angle `A` would imply im_target is the same shape as if im_reference was
rotated by `A`. 
Use `imrotate_function=imrotate_bin_<clockwise|counterclockwise>_<radians|degrees>` to get angles <clockwise|counterclockwise> in <radians|degrees>.
"""
function shape_difference_rotation(im_reference, im_target, test_angles; imrotate_function=imrotate_bin_clockwise_radians)
    imref_padded, imtarget_padded = pad_images(im_reference, im_target)
    shape_differences = Array{
        NamedTuple{(:angle, :shape_difference),Tuple{Float64,Float64}}
    }(
        undef, length(test_angles)
    )
    init_props = regionprops_table(label_components(im_reference))[1, :] # assumption only one object in image!
    idx = 1
    # r_init, c_init = compute_centroid(imref_padded, rounded=true)
    for angle in test_angles

        # rotate image back by angle
        imtarget_rotated = imrotate_function(imtarget_padded, -angle)

        im1, im2 = crop_to_shared_centroid(imref_padded, imtarget_rotated)

        # Check here that im1 and im2 sizes are the same
        if isequal(prod(size(im1)), prod(size(im2)))
            a_not_b = im1 .> 0 .&& isequal.(im2, 0)
            b_not_a = im2 .> 0 .&& isequal.(im1, 0)
            shape_difference = sum(a_not_b .|| b_not_a)
            shape_differences[idx] = (; angle, shape_difference)
        else
            @warn("Warning: shapes not equal\n")
            @warn(angle, size(im1), size(im2), "\n")
            shape_differences[idx] = (; angle, shape_difference=NaN)
        end
        idx += 1
    end
    return shape_differences
end


"""
Finds the image rotation angle in `test_angles` which minimizes the shape difference between `im_reference` and `im_target`.
The default test angles are evenly distributed in steps of π/36 rad (5º) around a full rotation,
ensuring that no angles are repeated (since -π rad == π rad),
and ordered so that smaller absolute angles which are positive will be returned in the event of a tie in the shape difference.
Use `imrotate_function=imrotate_bin_<clockwise|counterclockwise>_<radians|degrees>` to get angles <clockwise|counterclockwise> in <radians|degrees>.
"""

function register(
    im_reference,
    im_target;
    test_angles=sort(reverse(range(; start=-π, stop=π, step=π / 36)[1:(end-1)]); by=abs),
    imrotate_function=imrotate_bin_clockwise_radians,
)
    shape_differences = shape_difference_rotation(im_reference, im_target, test_angles; imrotate_function)
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    return best_match.angle
end


