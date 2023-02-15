# Helper functions
"""
    make_filename()

Makes default filename with timestamp.

"""
function make_filename()::String
    return "persisted_img-" * Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") * ".png"
end

function make_filename(fname::T, ext::T=".png")::T where T<:AbstractString
    return timestamp(fname) * ext
end

"""
    timestamp(fname)

Attach timestamp to `fname`.
"""
function timestamp(fname::String)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
    return fname * "-" * ts   
end

"""
    fname_ext_split(fname)

Split `"fname.ext"` into `"fname"` and `"ext"`.
"""
function fname_ext_split(fname::String)
    (name=fname[1:end-4],ext=fname[end-2:end])    
end

"""
    fname_ext_splice(fname, ext)

Join `"fname"` and `"ext"` with `'.'`.
"""
function fname_ext_splice(fname::String, ext::String)
    fname * '.' * ext
end


"""
    check_fname(fname)

Checks `fname` does not exist in current directory; throws an assertion if this condition is false.

# Arguments
- `fname`: String object or Symbol to a reference to a String representing a path.
"""
function check_fname(fname::Union{String,Symbol,Nothing}=nothing)::String
    if fname isa String # then use as filename
        check_name = fname
    elseif fname isa Symbol
        check_name = eval(fname) # get the object represented by the symbol
    elseif isnothing(fname) # nothing provided so make a filename
        check_name = make_filename()
    end

    # check name does not exist in wd
    isfile(check_name) && error("$check_name already exists in $(pwd())")
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
function bwdist(bwimg::AbstractArray{Bool})::AbstractArray{Float64}
    return Images.distance_transform(Images.feature_transform(bwimg))
end

"""
    padnhood(img, I, nhood)

Pad the matrix `img[nhood]` with zeros according to the position of `I` within the edges`img`.

Returns `img[nhood]` if `I` is not an edge index.
"""
function padnhood(img, I, nhood)
    # adaptive padding
    maxr, maxc = size(img)
    tofill = SizedMatrix{3,3}(zeros(Int, 3, 3))
    @views if I == CartesianIndex(1, 1) # top left corner`
        tofill[2:3, 2:3] = img[nhood]
    elseif I == CartesianIndex(maxr, 1) # bottom left corner 
        tofill[1:2, 2:3] = img[nhood]
    elseif I == CartesianIndex(1, maxc) # top right corner 
        tofill[2:3, 1:2] = img[nhood]
    elseif I == CartesianIndex(maxr, maxc) # bottom right corner 
        tofill[1:2, 1:2] = img[nhood]
    elseif I[1] == 1 # top edge (first row)
        tofill[2:3, 1:3] = img[nhood]
    elseif I[2] == 1 # left edge (first col)
        tofill[1:3, 2:3] = img[nhood]
    elseif I[1] == maxr # bottom edge (last row)
        tofill[1:2, 1:3] = img[nhood]
    elseif I[2] == maxc # right edge (last row)
        tofill[1:3, 1:2] = img[nhood]
    else
        tofill = img[nhood]
    end
    return tofill
end

"""
    _bin9todec(v)

Get decimal representation of a bit vector `v` with the leading bit at its leftmost posistion.

Example
```
julia> _bin9todec([0 0 0 0 0 0 0 0 0])    
0

julia> _bin9todec([1 1 1 1 1 1 1 1 1])    
511
```
"""
function _bin9todec(v::AbstractArray)::Int64
    return sum(vec(v) .* 2 .^ (0:(length(v) - 1)))
end

"""
    _operator_lut(I, img, nhood, lut1, lut2)

Look up the neighborhood `nhood` in lookup tables `lut1` and `lut2`.

Handles cases when the center of `nhood` is on the edge of `img` using data in `I`.
"""
function _operator_lut(
    I::CartesianIndex{2},
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
    lut1::Vector{Int64},
    lut2::Vector{Int64},
)::SVector{2, Int64}

    # corner pixels
    length(nhood) == 4 && return @SVector [false, 0]

    val = IceFloeTracker._bin9todec(_pad_handler(I, img, nhood)) + 1

    return @SVector [lut1[val], lut2[val]]
end

function _operator_lut(
    I::CartesianIndex{2},
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
    lut::Vector{T},
)::T where {T} # for bridge

    # corner pixels
    length(nhood) == 4 && return false # for bridge and some other operations like hbreak, branch

    return lut[_bin9todec(_pad_handler(I, img, nhood)) + 1]
end

function _pad_handler(I, img, nhood)
    (length(nhood) == 6) && return padnhood(img, I, nhood) # edge pixels
    return @view img[nhood]
end

"""
    check_2_tuple(pair::Tuple{Int64,Int64})

Check that the first element of `pair` is less than the second.
"""
function check_2_tuple(pair::Tuple{Int64,Int64})
    # assert the first area threshold is less than the second
    pair[1] >= pair[2] && throw(ArgumentError("The first element of $pair must be less than the second."))
    nothing
end

"""
    parse_2_tuple(s::String)

Parse a string of the form "a,b" into a tuple of Int64.
"""
function parse_2_tuple(s::String)::Tuple{Int64,Int64}
    # assert first and last characters are parentheses
    s[1] != '(' && throw(ArgumentError("The first character of `$s` must be a left parenthesis."))
    s[end] != ')' && throw(ArgumentError("The last character of `$s` must be a right parenthesis."))
    Tuple(parse.(Int64, split(s[2:end-1], ',')))
end

