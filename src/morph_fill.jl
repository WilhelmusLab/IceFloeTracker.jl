# Convert [1 1 1      [1 1 1
#          1 0 1  to   1 1 1
#          1 1 1]      1 1 1]

const __FILL_PATTERN__ = Bool[1 1 1; 1 0 1; 1 1 1]

function fill_operator(
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
)::Bool
    size(nhood) == (3, 3) && return img[nhood] == __FILL_PATTERN__
    return false
end

function morph_fill(img::T)::T where {T<:AbstractArray{Bool}}
    out = zeros(Bool, size(img))
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)

    for I in R
        if !img[I] # Check only zero pixels
            nhood = max(I_first, I - Δ):min(I_last, I + Δ)
            out[I] = fill_operator(img, nhood)
        end
    end

    return out .|| img
end
