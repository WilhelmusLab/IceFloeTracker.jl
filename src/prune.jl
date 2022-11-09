 """
    prune(img::AbstractArray{Bool})

Remove foreground pixels in the binary image `img` with degree less than 3 (fewer than 3 linear or diagonal neighboring `1`s).

# Examples
```jldoctest; setup = :(using IceFloeTracker)

julia> bw = zeros(Bool, 15, 15);

julia> bw[:,8] .= true; bw[8,:] .= true; bw[7:9, 7:9] .= true;

julia> [bw[i,i] = 1 for i=1:15]; # fill diagonal

julia> [bw[15-i,i+1] = 1 for i=0:14]; # fill the other diagonal

julia> bw
15×15 Matrix{Bool}:
 1  0  0  0  0  0  0  1  0  0  0  0  0  0  1
 0  1  0  0  0  0  0  1  0  0  0  0  0  1  0
 0  0  1  0  0  0  0  1  0  0  0  0  1  0  0
 0  0  0  1  0  0  0  1  0  0  0  1  0  0  0
 0  0  0  0  1  0  0  1  0  0  1  0  0  0  0
 0  0  0  0  0  1  0  1  0  1  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  0  0  0  0  0  0
 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
 0  0  0  0  0  0  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  1  0  1  0  1  0  0  0  0  0
 0  0  0  0  1  0  0  1  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  1  0  0  0  1  0  0  0
 0  0  1  0  0  0  0  1  0  0  0  0  1  0  0
 0  1  0  0  0  0  0  1  0  0  0  0  0  1  0
 1  0  0  0  0  0  0  1  0  0  0  0  0  0  1

julia> prune(bw)
15×15 Matrix{Bool}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  1  1  1  1  1  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
"""
function prune(img::T)::T where T<:AbstractArray{Bool}
    out = ones(Bool,size(img))
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)
    for I in R
        nhood = max(I_first, I-Δ):min(I_last, I+Δ)
        out[I] = isprunable(I, img, nhood)
    end
    return out
end

function isprunable(I::CartesianIndex{2},
                    img::AbstractArray{Bool}, nhood::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}})
    if img[I]
        return sum(img[nhood]) > 3
    else
        return false
    end
end
