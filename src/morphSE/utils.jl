
#=
    maybe_floattype(T)

To keep consistant with Base `diff` that Int array outputs Int array. It is sometimes
useful to only promote types for Bool and FixedPoint. In most of the time, `floattype`
should be the most reliable way.
=#
maybe_floattype(::Type{T}) where {T} = T
maybe_floattype(::Type{Bool}) = floattype(Bool)
# maybe_floattype(::Type{T}) where {T<:FixedPoint} = floattype(T)
# maybe_floattype(::Type{CT}) where {CT<:Color} = base_color_type(CT){maybe_floattype(eltype(CT))}

#=
A helper to eagerly check if particular function makes sense for `extreme_filter` semantics.
For instance, `max(::RGB, ::RGB)` is not well-defined and we should early throw errors so that
our users don't get encrypted error messages.
=#
require_select_function(f, ::Type{T}) where {T} = require_select_function(f, T, T)
function require_select_function(f, ::Type{T1}, ::Type{T2}) where {T1,T2}
    if !_is_select_function(f, T1, T2)
        hint = "does `f(x::T1, y::T2)` work as expected?"
        throw(
            ArgumentError(
                "function `$f` is not a well-defined select function on type `$T1` and `$T2`: $hint",
            ),
        )
    end
end
_is_select_function(f, ::Type{T}) where {T} = _is_select_function(f, T, T)

function _is_select_function(f, ::Type{T1}, ::Type{T2}) where {T1,T2}
    return _is_select_function_trial(f, T1, T2)
end
function _is_select_function(f, ::Type{T1}, ::Type{T2}) where {T1<:Real,T2<:Real}
    f in (min, max) && return true
    return _is_select_function_trial(f, T1, T2)
end
# function _is_select_function(f, ::Type{CT1}, ::Type{CT2}) where {CT1<:AbstractGray,CT2<:AbstractGray}
#     f in (min, max) && return true
#     return _is_select_function_trial(f, CT1, CT2)
# end
# function _is_select_function(f, ::Type{CT1}, ::Type{CT2}) where {CT1<:Colorant,CT2<:Colorant}
#     # min/max is not well-defined on generic color space
#     f in (min, max) && return false
#     return _is_select_function_trial(f, CT1, CT2)
# end
function _is_select_function_trial(f, ::Type{T1}, ::Type{T2}) where {T1,T2}
    # for generic case, just run a trial and see if it doesn't error
    try
        f(zero(T1), zero(T2))
        return true
    catch
        return false
    end
end

"""
    upper, lower = strel_split([T], se)

Split a symmetric structuring element into its upper and lower half parts based on its center point.

For each element `o` in `strel(CartesianIndex, upper)`, its negative `-o` is an element of `strel(CartesianIndex, lower)`. This function is not the inverse of [`strel_chain`](@ref).

The splited non-symmetric SE parts will be represented as array of `T`, where `T` is either a `Bool` or `CartesianIndex`. By default, `T = eltype(se)`.

```jldoctest strel_split; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> se = strel_diamond((3, 3))
3×3 SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> upper, lower = strel_split(se);

julia> upper
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  1  0
 1  1  0
 0  0  0

julia> lower
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  0  0
 0  1  1
 0  1  0
```

If the `se` is represented as displacement offset array, then the splitted result will also be displacement offset array:

```jldoctest strel_split
julia> se = strel(CartesianIndex, se)
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(0, -1)
 CartesianIndex(-1, 0)
 CartesianIndex(1, 0)
 CartesianIndex(0, 1)

julia> upper, lower = strel_split(se);

julia> upper
2-element Vector{CartesianIndex{2}}:
 CartesianIndex(0, -1)
 CartesianIndex(-1, 0)

julia> lower
2-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 0)
 CartesianIndex(0, 1)
```
"""
strel_split(se) = strel_split(eltype(se), se)
function strel_split(T, se)
    se = strel(Bool, se)
    require_symmetric_strel(se)
    R = LinearIndices(se)
    c = R[OffsetArrays.center(R)...]

    upper = copy(se)
    lower = copy(se)
    upper[(c + 1):end] .= false
    lower[begin:(c - 1)] .= false

    return strel(T, upper), strel(T, lower)
end

"""
    is_symmetric(se)

Check if a given structuring element array `se` is symmetric with respect to its center pixel.

More formally, this checks if `mask[I] == mask[-I]` for any valid `I ∈
CartesianIndices(mask)` in the connectivity mask representation `mask = strel(Bool, se)`.
"""
function is_symmetric(se::AbstractArray)
    # first check the axes, and then the values
    se = OffsetArrays.centered(strel(Bool, se))
    all(r -> first(r) == -last(r), axes(se)) || return false
    R = CartesianIndices(map(r -> 0:maximum(r), axes(se)))
    return all(R) do o
        @inbounds se[o] == se[-o]
    end
end
is_symmetric(se::SEBoxArray) = true
is_symmetric(se::SEDiamondArray) = true

#=
Some morphological operation only makes sense for symmetric structuring elements.
Here we provide a checker in spirit of Base.require_one_based_indexing.
=#
function require_symmetric_strel(se)
    return is_symmetric(se) || throw(
        ArgumentError("structuring element must be symmetric with respect to its center")
    )
end
