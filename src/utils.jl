# Helper functions
"""
    make_filename()

Makes default filename with timestamp.

"""
function make_filename()
    return "persisted_mask-" * Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") * ".png"
end

"""
    check_fname(fname)

Checks `fname` does not exist in current directory; throws an assertion if this condition is false.

# Arguments
- `fname`: String object or Symbol to a reference to a String representing a path.
"""
function check_fname(fname::Union{String,Symbol,Nothing}=nothing)
    if fname isa String # then use as filename
        check_name = fname
    elseif fname isa Symbol
        check_name = eval(fname) # get the object represented by the symbol
    elseif isnothing(fname) # nothing provided so make a filename
        check_name = make_filename()
    end

    # check name does not exist in wd
    @assert !isfile(check_name) "$check_name already exists in $(pwd())"
    return check_name
end

"""
    add_padding(img, style)

Extrapolate the image `img` according to the `style` specifications type. Returns the extrapolated image.

# Arguments
- `img`: Image to be padded.
- `style`: A supported type (such as `Pad` or `Fill`) representing the extrapolation style. See the relevant [documentation](https://juliaimages.org/latest/function_reference/#ImageFiltering) for details.

See also [`remove_padding`](@ref)
"""
function add_padding(img, style::Union{Pad,Fill})::Matrix
    return collect(Images.padarray(img, style))
end

"""
    remove_padding(paddedimg, border_spec)

Removes padding from the boundary of padded image `paddedimg` according to the border specification `border_spec` type. Returns the cropped image.

# Arguments
- `paddedimg`: Pre-padded image.
- `border_spec`: Type representing the style of padding (such as `Pad` or `Fill`) with which `paddedimg` is assumend to be pre-padded. Example: `Pad((1,2), (3,4))` specifies 1 row on the top, 2 columns on the left, 3 rows on the bottom, and 4 columns on the right boundary.

See also [`add_padding`](@ref)
"""
function remove_padding(paddedimg, border_spec::Union{Pad,Fill})::Matrix
    top, left = border_spec.lo
    bottom, right = border_spec.hi
    return paddedimg[(top + 1):(end - bottom), (left + 1):(end - right)]
end

"""
    imextendedmin(img)

Mimics MATLAB's imextendedmin function that computes the extended-minima transform, which is the regional minima of the H-minima transform. Regional minima are connected components of pixels with a constant intensity value. This function returns a transformed bitmatrix.

# Arguments
- `img`: image object
- `h`: suppress minima below this depth threshold
- `conn`: neighborhood connectivity; in 2D 1 = 4-neighborhood and 2 = 8-neighborhood
"""
function imextendedmin(img::AbstractArray; h::Int=2, conn::Int=2)::BitMatrix
    mask = ImageSegmentation.hmin_transform(img, h)
    mask_minima = Images.local_minima(mask; connectivity=conn)
    return Bool.(mask_minima)
end

"""
    bwdist(bwimg)

Distance transform for binary image `bwdist`.
"""
function bwdist(bwimg::AbstractArray{Bool})::T where {T<:AbstractArray{Float64}}
    return Images.distance_transform(
        Images.feature_transform(bwimg)
    )
end

"""
    padnhood(img, I, nhood)

Pad the matrix `img[nhood]` with zeros according to the position of `I` within the edges of `img`.

Returns `img[nhood]` if `I` is not an edge index.
"""
function padnhood(img::T, I::CartesianIndex{2},
    nhood::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}})::T where T<:AbstractArray{Bool}
        # adaptive padding
        maxr, maxc = size(img)
        tofill = zeros(Int,3,3);
        if I == CartesianIndex(1,1) # top left corner`
            tofill[2:3,2:3] = img[nhood]
        elseif I == CartesianIndex(maxr,1) # bottom left corner 
            tofill[1:2,2:3] = img[nhood]
        elseif I == CartesianIndex(1,maxc) # top right corner 
            tofill[2:3,1:2] = img[nhood]
        elseif I == CartesianIndex(maxr,maxc) # bottom right corner 
            tofill[1:2,1:2] = img[nhood]
        elseif I[1] == 1 # top edge (first row)
            tofill[2:3,1:3] = img[nhood]
        elseif I[2] == 1 # left edge (first col)
            tofill[1:3,2:3] = img[nhood]
        elseif I[1] == maxr # bottom edge (last row)
            tofill[1:2,1:3] = img[nhood]
        elseif I[2] == maxc # right edge (last row)
            tofill[1:3,1:2] = img[nhood]
        else
            tofill = img[nhood]
        end
        return tofill
    end
