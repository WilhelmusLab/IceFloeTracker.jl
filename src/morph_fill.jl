include("./lut/lutfill.jl")

const LUTFILL = make_lutfill()

function _fill_operator_lut(
    I::CartesianIndex{2},
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
)
    return _operator_lut(I, img, nhood, LUTFILL)
end

"""
    morph_fill(bw::T)::T where {T<:AbstractArray{Bool}}

Fill holes in binary image `bw` by setting 0-valued pixels to 1 if they are surrounded by 1-valued pixels.

# Examples

```jldoctest; setup = :(using IceFloeTracker)
julia> bw = Bool[
        0 0 0 0 0
        0 1 1 1 0
        0 1 0 1 0
        0 1 1 1 0
        0 0 0 0 0
    ];

julia> morph_fill(bw)
5×5 Matrix{Bool}:
 0  0  0  0  0
 0  1  1  1  0
 0  1  1  1  0
 0  1  1  1  0
 0  0  0  0  0
```
"""
function morph_fill(bw::T)::T where {T<:AbstractArray{Bool}}
    # TODO: see about implemting _filter using parallelization
    return _filter(bw, _fill_operator_lut)
end
