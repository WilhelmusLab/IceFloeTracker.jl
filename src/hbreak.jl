"""
    _make_hbreak_dict()

Build dict with the two versions of an H-connected 3x3 neighboorhood.

h1 =   [1 0 1
        1 1 1
        1 0 1]

h2 =   [1 1 1
        0 1 0
        1 1 1]
"""
function _make_hbreak_dict()::Dict{AbstractArray{Bool},Bool}
    h1 = trues(3, 3)
    h1[[1 3], 2] .= false
    h2 = trues(3, 3)
    h2[2, [1 3]] .= false
    d = Dict{AbstractArray{Bool},Bool}()
    d[h1], d[h2] = true, true
    return d
end

"""
    hbreak(img::AbstractArray{Bool})

Remove H-connected pixels in the binary image `img`. See also [`hbreak!`](@ref) for an inplace version of this function.

# Examples
```jldoctest; setup = :(using IceFloeTracker)

julia> h1 = trues(3,3); h1[[1 3], 2] .= false; h1
3×3 BitMatrix:
 1  0  1
 1  1  1
 1  0  1

julia> h2 = trues(3,3); h2[2, [1 3]] .= false; h2
3×3 BitMatrix:
 1  1  1
 0  1  0
 1  1  1

julia> hbreak!(h1); h1 # modify h1 inplace
3×3 BitMatrix:
 1  0  1
 1  0  1
 1  0  1

julia> hbreak(h2)
3×3 BitMatrix:
 1  1  1
 0  0  0
 1  1  1
```
"""
hbreak(img) = hbreak!(copy(img))

"""
    hbreak!(img::AbstractArray{Bool})

Inplace version of `hbreak`. See also [`hbreak`](@ref).
"""
function hbreak!(img::T)::T where {T<:AbstractArray{Bool}}
    hbreak_dict = _make_hbreak_dict()
    f = x -> get(hbreak_dict, x, false)
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)
    for I in R
        nhood = max(I_first, I - Δ):min(I_last, I + Δ)
        f(view(img, nhood)) ? img[I] = false : continue
    end
    return img
end
