greaterthan05(x) = x .> 0.5 # used for the image resize step and for binarizing images
function imrotate_bin(x, r)
    return greaterthan05(collect(imrotate(x, r, axes(x); method=BSpline(Constant()))))
end
function imrotate_bin_nocrop(x, r)
    return greaterthan05(collect(imrotate(x, r; method=BSpline(Constant()))))
end

# Functions used for the SD minimization
"""
Pad images by zeros based on the size of the larger of the two images.
"""
function pad_images(im1, im2)
    max1 = maximum(size(im1))
    max2 = maximum(size(im2))

    n = Int64(ceil(maximum([max1, max2])))
    im1_padded = collect(padarray(im1, Fill(0, (n, n), (n, n))))
    im2_padded = collect(padarray(im2, Fill(0, (n, n), (n, n))))
    return im1_padded, im2_padded
end

"""
Calculate the centroid of a binary image. If 'rounded', return the
nearest integer.
"""
function compute_centroid(im; rounded=false)
    xi = 0
    yi = 0
    R = sum(im .> 0)
    n, m = size(im)
    for ii in range(1, n)
        for jj in range(1, m)
            if im[ii, jj] > 0
                xi += ii
                yi += jj
            end
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
Computes the shape difference between im_reference and im_target for each angle (degrees) in test_angles.
The reference image is held constant, while the target image is rotated. The test_angles are interpreted
as the angle of rotation from target to reference, so to find the best match, we rotate the reverse
direction. A perfect match at angle A would imply im_target is the same shape as if im_reference was
rotated by A degrees. Use `mode=:counterclockwise` to get counterclockwise angles.
"""
function shape_difference_rotation(im_reference, im_target, test_angles; mode=:clockwise)
    imref_padded, imtarget_padded = pad_images(im_reference, im_target)
    shape_differences = Array{
        NamedTuple{(:angle, :shape_difference),Tuple{Float64,Float64}}
    }(
        undef, length(test_angles)
    )
    # shape_differences = zeros((length(test_angles), 2))
    init_props = regionprops_table(label_components(im_reference))[1, :] # assumption only one object in image!
    idx = 1
    # r_init, c_init = compute_centroid(imref_padded, rounded=true)
    for angle in test_angles

        if mode === :clockwise
            _angle = angle  # no-op
        elseif mode === :counterclockwise
            _angle = -angle  # image rotation algorithm works clockwise by default
        else
            throw("mode $(mode) not recognized")
        end

        # try rotating image back by angle
        imtarget_rotated = imrotate_bin(imtarget_padded, -_angle)

        im1, im2 = crop_to_shared_centroid(imref_padded, imtarget_rotated)

        # Check here that im1 and im2 sizes are the same
        # Could also add check that the images are nonempty
        # These checks could go inside the crop_to_shared_ccentroid function
        if isequal.(prod(size(im1)), prod(size(im2)))
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

function register(
    mask1,
    mask2;
    test_angles=sort(reverse(range(; start=-π, stop=π, step=π / 36)[1:(end-1)]); by=abs),
    mode=:clockwise,
)
    shape_differences = shape_difference_rotation(mask1, mask2, test_angles; mode)
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    return best_match.angle
end
