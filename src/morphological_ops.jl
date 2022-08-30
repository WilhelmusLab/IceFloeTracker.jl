"""
    bwperim(bwimg, conn::Int=4)

Locate the pixels at the boundary of objects in an binary image `bwimg` using connectivity `conn`.

# Arguments
- `bwimg`: Binary (black/white -- 1/0) image
- `conn`: connectivity. Either 4 or 8 for 4-pixel and 8-pixel connectivity, respectively.

# Examples

```jldoctest; setup = :(using IceFloeTracker)
julia> A = zeros(Int, 13, 16); A[2:6, 2:6] .= 1; A[4:8, 7:10] .= 1; A[10:12,13:15] .= 1; A[10:12,3:6] .= 1;

julia> A
13×16 Matrix{Int64}:
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
function bwperim(bwimg, conn::Int=4)
    # validate conn and define structuring element
    if typeof(conn) <: Int
        @assert conn in [4,8] "$conn is not a valid value for pixel connectivity. Permissible values are 4 (4-pixel connectivity) or 8 (8-pixel connectivity)."
        
        # Set se for the corresponding connectivity
        conn == 4 ?
            se = IceFloeTracker.strel_diamond((3,3)) : # 4 => diamond
            se = IceFloeTracker.strel_box((3,3))       # 8 => box
    end

    # work with BitArrays
    if typeof(bwimg) <: Matrix
        bwimg = BitArray(bwimg)
    end
    
    eroded_bwimg = IceFloeTracker.erode(bwimg, se)
    
    return bwimg .& .!eroded_bwimg
end
