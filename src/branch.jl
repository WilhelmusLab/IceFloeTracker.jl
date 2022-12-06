"""
    branch_candidates_func(nhood)

Filter `nhood` as candidate for branch point.
    
To be passed to the `make_lut` function.
"""
function branch_candidates_func(nhood::AbstractArray)::Bool
    nhood[2, 2] == 0 && return false
    sum(nhood) > 3 && return true
end

"""
    connected_background_count(nhood)

Second lut generator for neighbor transform with diamond strel (4-neighborhood).

To be passed to the `make_lut` function.

"""
function connected_background_count(nhood::AbstractArray)::Int64
    nhood[2, 2] != 0 && return maximum(label_components(.!Bool.(nhood)))
    return 0
end

"""
    make_lut(lutfunc::Function)

Generate lookup table (lut) for 3x3 neighborhoods according to `lutfunc`.
"""
function make_lut(lutfunc::Function)::Vector{Int}
    lut = vec(zeros(Int, 512))
    @inbounds @simd for i in 1:(2^9)
        v = parse.(Int, reverse(collect(bitstring(UInt16(i - 1))[8:end])))
        lut[i] = lutfunc(SMatrix{3,3}(v))
    end
    return lut
end

function _branch_operator_lut(
    I::CartesianIndex{2},
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
)
    lutbranchcandidates = make_lut(branch_candidates_func)
    lutbackcount4 = make_lut(connected_background_count)
    return _operator_lut(I, img, nhood, lutbranchcandidates, lutbackcount4)
end

"""
    _branch_filter(img::AbstractArray{Bool}, operator::Function)

Filter `img` with `operator`.
"""
function _branch_filter(
    img::T, operator::Function
)::Tuple{T,AbstractArray{Int64}} where {T<:AbstractArray{Bool}}
    C = zeros(Bool, size(img))
    B = zeros(Int, size(img))

    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)

    @inbounds @simd for I in R
        if img[I] # ones for branch
            nhood = max(I_first, I - Δ):min(I_last, I + Δ)
            C[I], B[I] = operator(I, img, nhood)
        end
    end
    return C, B
end

"""
    branch(img::AbstractArray{Bool})

Find branch points in skeletonized image `img` according to Definition 3 of [1].

[1] Arcelli, Carlo, and Gabriella Sanniti di Baja. "Skeletons of planar patterns." Machine Intelligence and Pattern Recognition. Vol. 19. North-Holland, 1996. 99-143.

"""
function branch(img::T)::T where {T<:AbstractArray{Bool}}
    # Get candidates C and 4-neighbor count
    C, B = _branch_filter(img, _branch_operator_lut)

    # Get non-end points E
    E = (B .!= 1)

    # Get final candidates FC (1st condition)
    FC = C .* E

    # Get pixels p with exactly 2 connected components (2nd condition)
    Vₚ = (B .== 2) .* E

    Vq = dilate((B .> 2) .* E)

    return FC .* .!(FC .* Vₚ .* Vq)
end
