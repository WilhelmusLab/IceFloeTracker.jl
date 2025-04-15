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
    r1, c1 = compute_centroid(im1; rounded=false)
    r2, c2 = compute_centroid(im2; rounded=false)

    n1, m1 = size(im1)
    n2, m2 = size(im2)
    new_halfn = minimum([r1, n1 - r1, r2, n2 - r2])
    new_halfm = minimum([c1, m1 - c1, c2, m2 - c2])

    _round(x) = round(Int32, x)
    im1_cropped = im1[
        _round(r1 - new_halfn):_round(r1 + new_halfn),
        _round(c1 - new_halfm):_round(c1 + new_halfm)
    ]
    im2_cropped = im2[
        _round(r2 - new_halfn):_round(r2 + new_halfn),
        _round(c2 - new_halfm):_round(c2 + new_halfm)
    ]

    return im1_cropped, im2_cropped
end

"""
Align images by padding so that the centroids of each image are on the edge of or within the same pixel.
"""
function align_centroids(im1::AbstractArray{Bool}, im2::AbstractArray{Bool})
    # Get the location of the pixel containing the centroids of im1 and im2 
    # in their current coordinate systems
    r1, c1 = Int64.(floor.(compute_centroid(im1; rounded=false)))
    r2, c2 = Int64.(floor.(compute_centroid(im2; rounded=false)))

    # Calculate the same centroid, but measured from the bottom right of each image
    s1, d1 = size(im1) .- (r1, c1) .+ 1
    s2, d2 = size(im2) .- (r2, c2) .+ 1

    # Calculate the new "common centroid" position in image coordinates
    rn, cn = (
        maximum([r1, r2]),
        maximum([c1, c2]),
    )
    # Calculate the new "reverse common centroid" position in image coordinates from the bottom right
    sn, dn = (
        maximum([s1, s2]),
        maximum([d1, d2]),
    )

    # For each image, we shift the pixel containing its centroid to the new centroid
    # by adding rn-ri rows padding at the top, and cn-ci columns at the left.
    # We ensure that the centroid is the same distance from the right border
    # by adding sn-si rows padding at the bottom and dn-di columns padding at the right
    # These need to be `collect`
    im1_padded = collect(padarray(im1, Fill(0, (rn - r1, cn - c1), (sn - s1, dn - d1))))
    im2_padded = collect(padarray(im2, Fill(0, (rn - r2, cn - c2), (sn - s2, dn - d2))))

    @assert floor.(compute_centroid(im1_padded; rounded=false)) == floor.(compute_centroid(im2_padded; rounded=false))

    return im1_padded, im2_padded
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
     

        # rotate image back by angle
        imtarget_rotated = imrotate_function(imtarget_padded, -angle)

        im1, im2 = align_centroids(imref_padded, imtarget_rotated)
        
        # Check here that im1 and im2 sizes are the same
        # should be guaranteed by "align_centroids"
        @assert size(im1) == size(im2)

        if isequal(size(im1), size(im2))
        a_not_b = im1 .> 0 .&& isequal.(im2, 0)
        b_not_a = im2 .> 0 .&& isequal.(im1, 0)
        shape_difference = sum(a_not_b .|| b_not_a)
        shape_differences[idx] = (; angle, shape_difference)
        idx += 1
    end
    return shape_differences
end


"""
The default registration angles are evenly distributed in steps of π/36 rad (5º) around a full rotation,
ensuring that no angles are repeated (since -π rad == π rad),
and ordered so that smaller absolute angles which are positive will be returned in the event of a tie in the shape difference.
"""
register_default_angles = sort(reverse(range(; start=-π, stop=π, step=π / 36)[1:(end-1)]); by=abs)

"""
Finds the image rotation angle in `test_angles` which minimizes the shape difference between `im_reference` and `im_target`.
The default test angles are shown in `register_default_angles`.
Use `imrotate_function=imrotate_bin_<clockwise|counterclockwise>_<radians|degrees>` to get angles <clockwise|counterclockwise> in <radians|degrees>.
"""
function register(
    im_reference,
    im_target;
    test_angles=register_default_angles,
    imrotate_function=imrotate_bin_clockwise_radians,
)
    shape_differences = shape_difference_rotation(im_reference, im_target, test_angles; imrotate_function)
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    return best_match.angle
end


