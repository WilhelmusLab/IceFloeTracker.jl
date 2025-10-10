import StaticArrays: @SVector, SVector, SizedMatrix

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
)::SVector{2,Int64}

    # corner pixels
    length(nhood) == 4 && return @SVector [false, 0]

    val = _bin9todec(_pad_handler(I, img, nhood)) + 1

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

function _pad_handler(I, img, nhood)
    (length(nhood) == 6) && return padnhood(img, I, nhood) # edge pixels
    return @view img[nhood]
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
