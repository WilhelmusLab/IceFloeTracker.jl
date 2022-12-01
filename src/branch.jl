"""
    branch_func(nhood)

Filter `nhood` as candidate for branch point. To be passed to the `make_lut` function
"""
function branch_func(nhood::AbstractArray)::Bool
    nhood[2,2] == 0 && return false
    sum(nhood) > 3 && return true
end

"""
    make_lut(lutfunc::Function)

Generate lookup table (lut) for 3x3 neighborhoods according to `lutfunc`.
"""
function make_lut(lutfunc::Function)::Vector{Bool}
    lut= vec(zeros(Bool,512))
        for i=1:2^9
            v = parse.(Int, reverse(collect(bitstring(UInt16(i-1))[8:end])))
            v = reshape(v, 3, 3)
            lut[i] = lutfunc(v)
        end
        return lut
end
    
"""
    connected_background_count(nhood::AbstractArray)::Int64

Second lut generator for neighbor transform with diamond strel (4-neighborhood).
"""
function connected_background_count(nhood::AbstractArray)::Int64
    nhood[2,2] != 0 && return maximum(label_components(.!Bool.(nhood)))
    return 0
end

lutbranch = make_lut(branch_func)


function _prelim_branch_operator_lut(I, img, nhood)
    return _operator_lut(I, img, nhood, lutbranch)
end

function _prelim_4conn_operator_lut(I, img, nhood)
    return _operator_lut(I, img, nhood, lutbc4in)
end

"""
    _filter(img::AbstractArray{Bool}, operator::Function)

Filter `img` with `operator`.
"""
function _filter(img::AbstractArray{Bool}, operator::Function) #::T where T<:AbstractArray{Bool}
    out = zeros(Int, size(img))
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)
    for I in R
        if  img[I] # ones for branch
            nhood = max(I_first, I-Δ):min(I_last, I+Δ)
            out[I]= operator(I, img, nhood)
        end
    end
    return out
end

