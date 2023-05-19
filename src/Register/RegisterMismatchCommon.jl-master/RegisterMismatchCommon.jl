module RegisterMismatchCommon

using ..RegisterCore

using ..CenterIndexedArrays, ImageCore

export correctbias!,
    nanpad,
    mismatch0,
    aperture_grid,
    allocate_mmarrays,
    default_aperture_width,
    truncatenoise!
export DimsLike,
    WidthLike,
    each_point,
    aperture_range,
    assertsamesize,
    tovec,
    mismatch,
    padsize,
    set_FFTPROD
export padranges, shiftrange, checksize_maxshift, register_translate

const DimsLike = Union{AbstractVector{Int},Dims}   # try to avoid these and just use Dims tuples for sizes
const WidthLike = Union{AbstractVector,Tuple}
FFTPROD = [2, 3]
set_FFTPROD(v) = global FFTPROD = v

function mismatch(
    fixed::AbstractArray{T},
    moving::AbstractArray{T},
    maxshift::DimsLike;
    normalization=:intensity,
) where {T<:AbstractFloat}
    return mismatch(T, fixed, moving, maxshift; normalization=normalization)
end
function mismatch(
    fixed::AbstractArray,
    moving::AbstractArray,
    maxshift::DimsLike;
    normalization=:intensity,
)
    return mismatch(Float32, fixed, moving, maxshift; normalization=normalization)
end

function mismatch_apertures(
    fixed::AbstractArray{T}, moving::AbstractArray{T}, args...; kwargs...
) where {T<:AbstractFloat}
    return mismatch_apertures(T, fixed, moving, args...; kwargs...)
end
function mismatch_apertures(fixed::AbstractArray, moving::AbstractArray, args...; kwargs...)
    return mismatch_apertures(Float32, fixed, moving, args...; kwargs...)
end

function mismatch_apertures(
    ::Type{T},
    fixed::AbstractArray,
    moving::AbstractArray,
    gridsize::DimsLike,
    maxshift::DimsLike;
    kwargs...,
) where {T}
    cs = coords_spatial(fixed)
    aperture_centers = aperture_grid(map(d -> size(fixed, d), cs), gridsize)
    aperture_width = default_aperture_width(fixed, gridsize)
    return mismatch_apertures(
        T, fixed, moving, aperture_centers, aperture_width, maxshift; kwargs...
    )
end

"""
`correctbias!(mm::MismatchArray)` replaces "suspect" mismatch
data with imputed data.  If each pixel in your camera has a different
bias, then matching that bias becomes an incentive to avoid
shifts.  Likewise, CMOS cameras tend to have correlated row/column
noise. These two factors combine to imply that `mm[i,j,...]` is unreliable
whenever `i` or `j` is zero.

Data are imputed by averaging the adjacent non-suspect values.  This
function works in-place, overwriting the original `mm`.
"""
function correctbias!(mm::MismatchArray{ND,N}, w=correctbias_weight(mm)) where {ND,N}
    T = eltype(ND)
    mxshift = maxshift(mm)
    Imax = CartesianIndex(mxshift)
    Imin = CartesianIndex(map(x -> -x, mxshift)::NTuple{N,Int})
    I1 = CartesianIndex(ntuple(d -> d > 2 ? 0 : 1, N)::NTuple{N,Int})  # only first 2 dims
    for I in eachindex(mm)
        if w[I] == 0
            mms = NumDenom{T}(0, 0)
            ws = zero(T)
            strt = max(Imin, I - I1)
            stop = min(Imax, I + I1)
            for J in CartesianIndices(ntuple(d -> strt[d]:stop[d], N))
                wJ = w[J]
                if wJ != 0
                    mms += wJ * mm[J]
                    ws += wJ
                end
            end
            mm[I] = mms / ws
        end
    end
    return mm
end

"`correctbias!(mms)` runs `correctbias!` on each element of an array-of-MismatchArrays."
function correctbias!(mms::AbstractArray{M}) where {M<:MismatchArray}
    for mm in mms
        correctbias!(mm)
    end
    return mms
end

function correctbias_weight(mm::MismatchArray{ND,N}) where {ND,N}
    T = eltype(ND)
    w = CenterIndexedArray(ones(T, size(mm)))
    for I in eachindex(mm)
        anyzero = false
        for d in 1:min(N, 2)   # only first 2 dims
            anyzero |= I[d] == 0
        end
        if anyzero
            w[I] = 0
        end
    end
    return w
end

"""
`fixedpad, movingpad = nanpad(fixed, moving)` will pad `fixed` and/or
`moving` with NaN as needed to ensure that `fixedpad` and `movingpad`
have the same size.
"""
function nanpad(fixed, moving)
    ndims(fixed) == ndims(moving) ||
        error("fixed and moving must have the same dimensionality")
    if size(fixed) == size(moving)
        return fixed, moving
    end
    rng = map(d -> 1:max(size(fixed, d), size(moving, d)), 1:ndims(fixed))
    T = promote_type(eltype(fixed), eltype(moving))
    return get(fixed, rng, nanval(T)), get(moving, rng, nanval(T))
end

nanval(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
nanval(::Type{T}) where {T} = convert(Float32, NaN)

"""
`mm0 = mismatch0(fixed, moving, [normalization])` computes the
"as-is" mismatch between `fixed` and `moving`, without any shift.
`normalization` may be either `:intensity` (the default) or `:pixels`.
"""
function mismatch0(
    fixed::AbstractArray{Tf,N}, moving::AbstractArray{Tm,N}; normalization=:intensity
) where {Tf,Tm,N}
    size(fixed) == size(moving) || throw(
        DimensionMismatch(
            "Size $(size(fixed)) of fixed is not equal to size $(size(moving)) of moving",
        ),
    )
    return _mismatch0(
        zero(Float64), zero(Float64), fixed, moving; normalization=normalization
    )
end

function _mismatch0(
    num::T,
    denom::T,
    fixed::AbstractArray{Tf,N},
    moving::AbstractArray{Tm,N};
    normalization=:intensity,
) where {T,Tf,Tm,N}
    if normalization == :intensity
        for i in eachindex(fixed, moving)
            vf = T(fixed[i])
            vm = T(moving[i])
            if isfinite(vf) && isfinite(vm)
                num += (vf - vm)^2
                denom += vf^2 + vm^2
            end
        end
    elseif normalization == :pixels
        for i in eachindex(fixed, moving)
            vf = T(fixed[i])
            vm = T(moving[i])
            if isfinite(vf) && isfinite(vm)
                num += (vf - vm)^2
                denom += 1
            end
        end
    else
        error("Normalization $normalization not recognized")
    end
    return NumDenom(num, denom)
end

"""
`mm0 = mismatch0(mms)` computes the "as-is"
mismatch between `fixed` and `moving`, without any shift.  The
mismatch is represented in `mms` as an aperture-wise
Arrays-of-MismatchArrays.
"""
function mismatch0(mms::AbstractArray{M}) where {M<:MismatchArray}
    mm0 = eltype(M)(0, 0)
    cr = eachindex(first(mms))
    z = first(cr) + last(cr)  # all-zeros CartesianIndex
    for mm in mms
        mm0 += mm[z]
    end
    return mm0
end

"""
`ag = aperture_grid(ssize, gridsize)` constructs a uniformly-spaced
grid of aperture centers.  The grid has size `gridsize`, and is
constructed for an image of spatial size `ssize`.  Along each
dimension the first and last elements are at the image corners.
"""
function aperture_grid(ssize::Dims{N}, gridsize) where {N}
    if length(gridsize) != N
        if length(gridsize) == N - 1
            @info(
                "ssize and gridsize disagree; possible fix is to use a :time axis (AxisArrays) for the image"
            )
        end
        error("ssize and gridsize must have the same length, got $ssize and $gridsize")
    end
    grid = Array{NTuple{N,Float64},N}(undef, (gridsize...,))
    centers = map(
        i -> if gridsize[i] > 1
            collect(range(1; stop=ssize[i], length=gridsize[i]))
        else
            [(ssize[i] + 1) / 2]
        end,
        1:N,
    )
    for I in CartesianIndices(size(grid))
        grid[I] = ntuple(i -> centers[i][I[i]], N)
    end
    return grid
end

"""
`mms = allocate_mmarrays(T, gridsize, maxshift)` allocates storage for
aperture-wise mismatch computation. `mms` will be an
Array-of-MismatchArrays with element type `NumDenom{T}` and half-size
`maxshift`. `mms` will be an array of size `gridsize`. This syntax is
recommended when your apertures are centered at points of a grid.

`mms = allocate_mmarrays(T, aperture_centers, maxshift)` returns `mms`
with a shape that matches that of `aperture_centers`. The centers can
in general be provided as an vector-of-tuples, vector-of-vectors, or a
matrix with each point in a column.  If your centers are arranged in a
rectangular grid, you can use an `N`-dimensional array-of-tuples (or
array-of-vectors) or an `N+1`-dimensional array with the center
positions specified along the first dimension.  (But you may find the
`gridsize` syntax to be simpler.)
"""
function allocate_mmarrays(
    ::Type{T}, aperture_centers::AbstractArray{C}, maxshift
) where {T,C<:Union{AbstractVector,Tuple}}
    isempty(aperture_centers) && error("aperture_centers is empty")
    N = length(first(aperture_centers))
    sz = map(x -> 2 * x + 1, maxshift)
    mm = MismatchArray(T, sz...)
    mms = Array{typeof(mm)}(undef, size(aperture_centers))
    f = true
    for i in eachindex(mms)
        if f
            mms[i] = mm
            f = false
        else
            mms[i] = MismatchArray(T, sz...)
        end
    end
    return mms
end

function allocate_mmarrays(
    ::Type{T}, aperture_centers::AbstractArray{R}, maxshift
) where {T,R<:Real}
    N = ndims(aperture_centers) - 1
    mms = Array{MismatchArray{T,N}}(undef, size(aperture_centers)[2:end])
    sz = map(x -> 2 * x + 1, maxshift)
    for i in eachindex(mms)
        mms[i] = MismatchArray(T, sz...)
    end
    return mms
end

function allocate_mmarrays(::Type{T}, gridsize::NTuple{N,Int}, maxshift) where {T<:Real,N}
    mms = Array{MismatchArray{NumDenom{T},N}}(undef, gridsize)
    sz = map(x -> 2 * x + 1, maxshift)
    for i in eachindex(mms)
        mms[i] = MismatchArray(T, sz...)
    end
    return mms
end

struct ContainerIterator{C}
    data::C
end

Base.iterate(iter::ContainerIterator) = iterate(iter.data)
Base.iterate(iter::ContainerIterator, state) = iterate(iter.data, state)

struct FirstDimIterator{A<:AbstractArray,R<:CartesianIndices}
    data::A
    rng::R

    function FirstDimIterator{A,R}(data::A) where {A,R}
        return new{A,R}(data, CartesianIndices(Base.tail(size(data))))
    end
end
function FirstDimIterator(A::AbstractArray)
    return FirstDimIterator{typeof(A),typeof(CartesianIndices(Base.tail(size(A))))}(A)
end
function Base.iterate(iter::FirstDimIterator)
    isempty(iter.rng) && return nothing
    index, state = iterate(iter.rng)
    return iter.data[:, index], state
end
function Base.iterate(iter::FirstDimIterator, state)
    state == last(iter.rng) && return nothing
    index, state = iterate(iter.rng, state)
    return iter.data[:, index], state
end

"""
`iter = each_point(points)` yields an iterator `iter` over all the
points in `points`. `points` may be represented as an
AbstractArray-of-tuples or -AbstractVectors, or may be an
`AbstractArray` where each point is represented along the first
dimension (e.g., columns of a matrix).
"""
function each_point(
    aperture_centers::AbstractArray{C}
) where {C<:Union{AbstractVector,Tuple}}
    return ContainerIterator(aperture_centers)
end

function each_point(aperture_centers::AbstractArray{R}) where {R<:Real}
    return FirstDimIterator(aperture_centers)
end

"""
`rng = aperture_range(center, width)` returns a tuple of
`UnitRange{Int}`s that, for dimension `d`, is centered on `center[d]`
and has width `width[d]`.
"""
function aperture_range(center, width)
    length(center) == length(width) || error("center and width must have the same length")
    return ntuple(
        d -> leftedge(center[d], width[d]):rightedge(center[d], width[d]), length(center)
    )
end

"""
`aperturesize = default_aperture_width(img, gridsize, [overlap])`
calculates the aperture width for a regularly-spaced grid of aperture
centers with size `gridsize`.  Apertures that are adjacent along
dimension `d` may overlap by a number pixels specified by
`overlap[d]`; the default value is 0.  For non-negative `overlap`, the
collection of apertures will yield full coverage of the image.
"""
function default_aperture_width(
    img, gridsize::DimsLike, overlap::DimsLike=zeros(Int, sdims(img))
)
    sc = coords_spatial(img)
    length(sc) == length(gridsize) == length(overlap) || error(
        "gridsize and overlap must have length equal to the number of spatial dimensions in img",
    )
    for i in 1:length(sc)
        if gridsize[i] > size(img, sc[i])
            error(
                "gridsize $gridsize is too large, given the size $(size(img)[sc]) of the image",
            )
        end
    end
    gsz1 = max.(1, [gridsize...] .- 1)
    gflag = [gridsize...] .> 1
    return tuple(
        (([map(d -> size(img, d), sc)...] - gflag) ./ gsz1 + 2 * [overlap...] .* gflag)...
    )
end

"""
`truncatenoise!(mm, thresh)` zeros out any entries of the
MismatchArray `mm` whose `denom` values are less than `thresh`.
"""
function truncatenoise!(mm::AbstractArray{NumDenom{T}}, thresh::Real) where {T<:Real}
    for I in eachindex(mm)
        if mm[I].denom <= thresh
            mm[I] = NumDenom{T}(0, 0)
        end
    end
    return mm
end

function truncatenoise!(mms::AbstractArray{A}, thresh::Real) where {A<:MismatchArray}
    for i in 1:length(denoms)
        truncatenoise!(mms[i], thresh)
    end
    return nothing
end

"""
`shift = register_translate(fixed, moving, maxshift, [thresh])`
computes the integer-valued translation which best aligns images
`fixed` and `moving`. All shifts up to size `maxshift` are considered.
Optionally specify `thresh`, the fraction (0<=thresh<=1) of overlap
required between `fixed` and `moving` (default 0.25).
"""
function register_translate(fixed, moving, maxshift, thresh=nothing)
    mm = mismatch(fixed, moving, maxshift)
    _, denom = separate(mm)
    if thresh == nothing
        thresh = 0.25maximum(denom)
    end
    return indmin_mismatch(mm, thresh)
end

function checksize_maxshift(A::AbstractArray, maxshift)
    ndims(A) == length(maxshift) || error(
        "Array is $(ndims(A))-dimensional, but maxshift has length $(length(maxshift))"
    )
    for i in 1:ndims(A)
        size(A, i) == 2 * maxshift[i] + 1 || error(
            "Along dimension $i, the output size $(size(A,i)) does not agree with maxshift[$i] = $(maxshift[i])",
        )
    end
    return nothing
end

function padranges(blocksize, maxshift)
    padright = [maxshift...]
    transformdims = findall(padright .> 0)
    paddedsz = [blocksize...] + 2 * padright
    for i in transformdims
        # Pick a size for which one can efficiently calculate ffts
        padright[i] += padsize(blocksize, maxshift, i) - paddedsz[i]
    end
    return rng = UnitRange{Int}[
        (1 - maxshift[i]):(blocksize[i] + padright[i]) for i in 1:length(blocksize)
    ]
end

function padsize(blocksize::Dims{N}, maxshift::Dims{N}) where {N}
    return map(padsize, blocksize, maxshift, ntuple(identity, Val(N)))
end

function padsize(blocksize::Int, maxshift::Int, dim::Int)
    p = blocksize + 2maxshift
    return maxshift > 0 ? (dim == 1 ? nextpow(2, p) : nextprod(FFTPROD, p)) : p   # we won't FFT along dimensions with maxshift 0
end

function assertsamesize(A, B)
    if !issamesize(A, B)
        error("Arrays are not the same size")
    end
end

function issamesize(A::AbstractArray, B::AbstractArray)
    n = ndims(A)
    ndims(B) == n || return false
    for i in 1:n
        axes(A, i) == axes(B, i) || return false
    end
    return true
end

function issamesize(A::AbstractArray, indices)
    n = ndims(A)
    length(indices) == n || return false
    for i in 1:n
        size(A, i) == length(indices[i]) || return false
    end
    return true
end

# This yields the _effective_ overlap, i.e., sets to zero if gridsize==1 along a coordinate
# imgssz = image spatial size
function computeoverlap(imgssz, blocksize, gridsize)
    gsz1 = max(1, [gridsize...] .- 1)
    tmp = [imgssz...] ./ gsz1
    return blocksize - [ceil(Int, x) for x in tmp]
end

leftedge(center, width) = ceil(Int, center - width / 2)
rightedge(center, width) = leftedge(center + width, width) - 1

# These avoid making a copy if it's not necessary
tovec(v::AbstractVector) = v
tovec(v::Tuple) = [v...]

shiftrange(r, s) = r .+ s

### Utilities for unsafe indexing of views
# TODO: redesign this whole thing to be safer?
using Base: ViewIndex, to_indices, unsafe_length, index_shape, tail

@inline function extraunsafe_view(V::SubArray{T,N}, I::Vararg{ViewIndex,N}) where {T,N}
    idxs = unsafe_reindex(V, V.indices, to_indices(V, I))
    return SubArray(V.parent, idxs)
end

function get_index_wo_boundcheck(r::AbstractUnitRange, s::AbstractUnitRange{<:Integer}) # BlockRegistration issue #36
    f = first(r)
    st = oftype(f, f + first(s) - 1)
    return range(st; length=length(s))
end

function unsafe_reindex(
    V, idxs::Tuple{UnitRange,Vararg{Any}}, subidxs::Tuple{UnitRange,Vararg{Any}}
)
    (Base.@_propagate_inbounds_meta;
    @inbounds new1 = get_index_wo_boundcheck(idxs[1], subidxs[1]);
    (new1, unsafe_reindex(V, tail(idxs), tail(subidxs))...))
end

unsafe_reindex(V, idxs, subidxs) = Base.reindex(V, idxs, subidxs)

### Deprecations

function padsize(blocksize, maxshift)
    Base.depwarn(
        "padsize(::$(typeof(blocksize)), ::$(typeof(maxshift)) is deprecated, use Dims-tuples instead",
        :padsize,
    )
    sz = Vector{Int}(undef, length(blocksize))
    return padsize!(sz, blocksize, maxshift)
end

function padsize!(sz::Vector, blocksize, maxshift)
    n = length(blocksize)
    for i in 1:n
        sz[i] = padsize(blocksize, maxshift, i)
    end
    return sz
end

function padsize(blocksize, maxshift, dim)
    m = maxshift[dim]
    p = blocksize[dim] + 2m
    return m > 0 ? (dim == 1 ? nextpow(2, p) : nextprod(FFTPROD, p)) : p   # we won't FFT along dimensions with maxshift[i]==0
end

end #module
