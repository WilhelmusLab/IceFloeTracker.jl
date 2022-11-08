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
function _bin9todec(v)
    return sum(vec(v) .* 2 .^ (0:length(v)-1))
end

function operator_lut(I, img, nhood, lut) # for bridge
    # corner pixels
    if length(nhood) == 4
        return false
    elseif length(nhood) == 6
        filled = padnhood(img, I, nhood)
    else
        filled = img[nhood]
    end
    return lut[_bin9todec(filled)+1]
end

function bridge_operator_lut(I, img, nhood)
    lutbridge = vec(readdlm("src/lut/lutbridge.csv", ',', Bool))
    return operator_lut(I, img, nhood, lutbridge)
end

function bridge_filter(img, operator)
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

function bridge(bw)
    if eltype(bw) != Bool
        bw = Bool.(bw)
    end
    bridge_filter(bw, bridge_operator_lut)
end
