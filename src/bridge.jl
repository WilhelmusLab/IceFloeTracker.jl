include("./lut/lutbridge.jl")

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
    return sum(vec(v) .* 2 .^ (0:length(v)-1))
end

function _operator_lut(I::CartesianIndex{2}, img::AbstractArray{Bool},
    nhood::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}},
    lut::Vector{T})::T where T # for bridge
    # corner pixels
    if length(nhood) == 4
        return false # false for bridge and some other operations like hbreak, branch
    elseif length(nhood) == 6 # edge pixels
        filled = padnhood(img, I, nhood)
    else # interior pixels
        filled = img[nhood]
    end
    return lut[_bin9todec(filled)+1]
end

function _bridge_operator_lut(I::CartesianIndex{2}, img::AbstractArray{Bool},
    nhood::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}})
    lutbridge = make_lutbridge()
    return _operator_lut(I, img, nhood, lutbridge)
end

function _bridge_filter(img::T, operator::Function)::T where T<:AbstractArray{Bool}
    out = zeros(Bool,size(img))
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)
    for I in R
        if !img[I] # zero pixels only
            nhood = max(I_first, I-Δ):min(I_last, I+Δ)        
            out[I]= operator(I, img, nhood)
        end
    end
    return out .|| img
end

"""
    bridge(bw)

Set 0-valued pixels to 1 if they have two nonzero neighbors that are not connected. Note the following exceptions:

0 0 0           0 0 0
1 0 1  becomes  1 1 1
0 0 0           0 0 0

1 0 1           1 1 1
0 0 0  becomes  0 0 0
0 0 0           0 0 0

The same applies to all their corresponding rotations.

# Examples

```jldoctest; setup = :(using IceFloeTracker)
julia> bw = [0 0 0; 0 0 0; 1 0 1]
3×3 Matrix{Int64}:
 0  0  0
 0  0  0
 1  0  1

julia> bridge(bw)
3×3 BitMatrix:
 0  0  0
 0  0  0
 1  1  1

julia> bw = [1 0 0; 1 0 1; 0 0 1]
3×3 Matrix{Int64}:
 1  0  0
 1  0  1
 0  0  1

julia> bridge(bw)
3×3 BitMatrix:
 1  1  0
 1  1  1
 0  1  1

 julia> bw = [1 0 1; 0 0 0; 1 0 1]
3×3 Matrix{Int64}:
 1  0  1
 0  0  0
 1  0  1

julia> bridge(bw)
3×3 BitMatrix:
 1  1  1
 1  1  1
 1  1  1
```
"""
function bridge(bw::T)::T where T<:AbstractArray{Bool}
    _bridge_filter(bw, _bridge_operator_lut)
end
