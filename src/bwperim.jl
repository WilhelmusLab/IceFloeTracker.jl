"""
    bwperim(bwimg)

Locate the pixels at the boundary of objects in an binary image `bwimg` using 8-pixel connectivity.

# Arguments
- `bwimg`: Binary (black/white -- 1/0) image

# Examples

```jldoctest; setup = :(using IceFloeTracker)
julia> A = zeros(Bool, 13, 16); A[2:6, 2:6] .= 1; A[4:8, 7:10] .= 1; A[10:12,13:15] .= 1; A[10:12,3:6] .= 1;

julia> A
13×16 Matrix{Bool}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

 julia> bwperim(A)
13×16 Matrix{Bool}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
 0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0
 0  1  0  0  0  0  1  1  1  1  0  0  0  0  0  0
 0  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0
 0  1  1  1  1  1  0  0  0  1  0  0  0  0  0  0
 0  0  0  0  0  0  1  0  0  1  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
```
"""
function bwperim(bwimg::T) where {T<:AbstractMatrix{Bool}}
    eroded_bwimg = erode(bwimg)
    return bwimg .& .!eroded_bwimg
end
