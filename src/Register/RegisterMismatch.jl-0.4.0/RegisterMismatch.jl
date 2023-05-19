module RegisterMismatch

import Base: copy, eltype, isnan, ndims

using ImageCore
using RFFT, FFTW
using ..RegisterCore, PaddedViews, MappedArrays
using Printf
using Reexport
@reexport using ..RegisterMismatchCommon
import ..RegisterMismatchCommon: mismatch0, mismatch, mismatch_apertures

export
    CMStorage,
    fillfixed!,
    mismatch0,
    mismatch,
    mismatch!,
    mismatch_apertures,
    mismatch_apertures!

"""
The major types and functions exported are:

- `mismatch` and `mismatch!`: compute the mismatch between two images
- `mismatch_apertures` and `mismatch_apertures!`: compute apertured mismatch between two images
- `mismatch0`: simple direct mismatch calculation with no shift
- `nanpad`: pad the smaller image with NaNs
- `highpass`: highpass filter an image
- `correctbias!`: replace corrupted mismatch data (due to camera bias inhomogeneity) with imputed data
- `truncatenoise!`: threshold mismatch computation to prevent problems from roundoff
- `aperture_grid`: create a regular grid of apertures
- `allocate_mmarrays`: create storage for output of `mismatch_apertures!`
- `CMStorage`: a type that facilitates re-use of intermediate storage during registration computations
"""
RegisterMismatch

FFTW.set_num_threads(min(Sys.CPU_THREADS, 8))
set_FFTPROD([2, 3])

mutable struct NanCorrFFTs{T<:AbstractFloat,N,RCType<:RCpair{T,N}}
    I0::RCType
    I1::RCType
    I2::RCType
end

copy(x::NanCorrFFTs) = NanCorrFFTs(copy(x.I0), copy(x.I1), copy(x.I2))

"""
    CMStorage{T}(undef, aperture_width, maxshift; flags=FFTW.ESTIMATE, timelimit=Inf, display=true)

Prepare for FFT-based mismatch computations over domains of size `aperture_width`, computing the
mismatch up to shifts of size `maxshift`.  The keyword arguments allow you to control the planning
process for the FFTs.
"""
mutable struct CMStorage{T<:AbstractFloat,N,RCType<:RCpair{T,N},FFT<:Function,IFFT<:Function}
    aperture_width::Vector{Float64}
    maxshift::Vector{Int}
    getindices::Vector{UnitRange{Int}}   # indices for pulling padded data, in source-coordinates
    padded::Array{T,N}
    fixed::NanCorrFFTs{T,N,RCType}
    moving::NanCorrFFTs{T,N,RCType}
    buf1::RCType
    buf2::RCType
    # the next two store the result of calling plan_fft! and plan_ifft!
    fftfunc!::FFT
    ifftfunc!::IFFT
    shiftindices::Vector{Vector{Int}} # indices for performing fftshift & snipping from -maxshift:maxshift
end

function CMStorage{T,N}(::UndefInitializer, aperture_width::NTuple{N,<:Real}, maxshift::Dims{N}; flags=FFTW.ESTIMATE, timelimit=Inf, display=true) where {T,N}
    blocksize = map(x -> ceil(Int, x), aperture_width)
    padsz = padsize(blocksize, maxshift)
    padded = Array{T}(undef, padsz)
    getindices = padranges(blocksize, maxshift)
    maxshiftv = [maxshift...]
    region = findall(maxshiftv .> 0)
    fixed = NanCorrFFTs(RCpair{T}(undef, padsz, region), RCpair{T}(undef, padsz, region), RCpair{T}(undef, padsz, region))
    moving = NanCorrFFTs(RCpair{T}(undef, padsz, region), RCpair{T}(undef, padsz, region), RCpair{T}(undef, padsz, region))
    buf1 = RCpair{T}(undef, padsz, region)
    buf2 = RCpair{T}(undef, padsz, region)
    tcalib = 0
    if display && flags != FFTW.ESTIMATE
        print("Planning FFTs (maximum $timelimit seconds)...")
        flush(stdout)
        tcalib = time()
    end
    fftfunc = plan_rfft!(fixed.I0, flags=flags, timelimit=timelimit / 2)
    ifftfunc = plan_irfft!(fixed.I0, flags=flags, timelimit=timelimit / 2)
    if display && flags != FFTW.ESTIMATE
        dt = time() - tcalib
        @printf("done (%.2f seconds)\n", dt)
    end
    shiftindices = Vector{Int}[[size(padded, i) .+ (-maxshift[i]+1:0); 1:maxshift[i]+1] for i = 1:length(maxshift)]
    CMStorage{T,N,typeof(buf1),typeof(fftfunc),typeof(ifftfunc)}(Float64[aperture_width...], maxshiftv, getindices, padded, fixed, moving, buf1, buf2, fftfunc, ifftfunc, shiftindices)
end

CMStorage{T}(::UndefInitializer, aperture_width::NTuple{N,<:Real}, maxshift::Dims{N}; kwargs...) where {T<:Real,N} =
    CMStorage{T,N}(undef, aperture_width, maxshift; kwargs...)

eltype(cms::CMStorage{T,N}) where {T,N} = T
ndims(cms::CMStorage{T,N}) where {T,N} = N

"""
    mm = mismatch([T], fixed, moving, maxshift; normalization=:intensity)

Compute the mismatch between `fixed` and
`moving` as a function of translations (shifts) up to size `maxshift`.
Optionally specify the element-type of the mismatch arrays (default
`Float32` for Integer- or FixedPoint-valued images) and the
normalization scheme (`:intensity` or `:pixels`).

`fixed` and `moving` must have the same size; you can pad with
`NaN`s as needed. See `nanpad`.
"""
function mismatch(::Type{T}, fixed::AbstractArray, moving::AbstractArray, maxshift::DimsLike; normalization=:intensity) where {T<:Real}
    msz = 2 .* maxshift .+ 1
    mm = MismatchArray(T, msz...)
    cms = CMStorage{T}(undef, size(fixed), maxshift)
    fillfixed!(cms, fixed)
    erng = shiftrange.((cms.getindices...,), first.(axes(fixed)) .- 1)  # expanded rng
    mpad = PaddedView(convert(T, NaN), of_eltype(T, moving), erng)
    mismatch!(mm, cms, mpad, normalization=normalization)
    return mm
end

"""
`mismatch!(mm, cms, moving; [normalization=:intensity])`
computes the mismatch as a function of shift, storing the result in
`mm`. The `fixed` image has been prepared in `cms`, a `CMStorage` object.
"""
function mismatch!(mm::MismatchArray, cms::CMStorage, moving::AbstractArray; normalization=:intensity)
    # Pad the moving snippet using any available data, including
    # regions that might be in the parent Array but are not present
    # within the boundaries of the SubArray. Use NaN only for pixels
    # truly lacking data.
    checksize_maxshift(mm, cms.maxshift)
    copyto!(cms.padded, CartesianIndices(cms.padded), moving, CartesianIndices(moving))
    fftnan!(cms.moving, cms.padded, cms.fftfunc!)
    # Compute the mismatch
    f0 = complex(cms.fixed.I0)
    f1 = complex(cms.fixed.I1)
    f2 = complex(cms.fixed.I2)
    m0 = complex(cms.moving.I0)
    m1 = complex(cms.moving.I1)
    m2 = complex(cms.moving.I2)
    tnum = complex(cms.buf1)
    tdenom = complex(cms.buf2)
    if normalization == :intensity
        @inbounds Threads.@threads for i = 1:length(tnum)
            c = 2 * conj(f1[i]) * m1[i]
            q = conj(f2[i]) * m0[i] + conj(f0[i]) * m2[i]
            tdenom[i] = q
            tnum[i] = q - c
        end
    elseif normalization == :pixels
        @inbounds Threads.@threads for i = 1:length(tnum)
            f0i, m0i = f0[i], m0[i]
            tdenom[i] = conj(f0i) * m0i
            tnum[i] = conj(f2[i]) * m0i - 2 * conj(f1[i]) * m1[i] + conj(f0i) * m2[i]
        end
    else
        error("normalization $normalization not recognized")
    end
    cms.ifftfunc!(cms.buf1)
    cms.ifftfunc!(cms.buf2)
    copyto!(mm, (view(real(cms.buf1), cms.shiftindices...), view(real(cms.buf2), cms.shiftindices...)))
    mm
end

"""
`mms = mismatch_apertures([T], fixed, moving, gridsize, maxshift;
[normalization=:pixels], [flags=FFTW.MEASURE], kwargs...)` computes
the mismatch between `fixed` and `moving` over a regularly-spaced grid
of aperture centers, effectively breaking the images up into
chunks. The maximum-allowed shift in any aperture is `maxshift`.

`mms = mismatch_apertures([T], fixed, moving, aperture_centers,
aperture_width, maxshift; kwargs...)` computes the mismatch between
`fixed` and `moving` over a list of apertures of size `aperture_width`
at positions defined by `aperture_centers`.

`fixed` and `moving` must have the same size; you can pad with `NaN`s
as needed to ensure this.  You can optionally specify the real-valued
element type mm; it defaults to the element type of `fixed` and
`moving` or, for Integer- or FixedPoint-valued images, `Float32`.

On output, `mms` will be an Array-of-MismatchArrays, with the outer
array having the same "grid" shape as `aperture_centers`.  The centers
can in general be provided as an vector-of-tuples, vector-of-vectors,
or a matrix with each point in a column.  If your centers are arranged
in a rectangular grid, you can use an `N`-dimensional array-of-tuples
(or array-of-vectors) or an `N+1`-dimensional array with the center
positions specified along the first dimension. See `aperture_grid`.
"""
function mismatch_apertures(::Type{T},
    fixed::AbstractArray,
    moving::AbstractArray,
    aperture_centers::AbstractArray,
    aperture_width::WidthLike,
    maxshift::DimsLike;
    normalization=:pixels,
    flags=FFTW.MEASURE,
    kwargs...) where {T}
    nd = sdims(fixed)
    (length(aperture_width) == nd && length(maxshift) == nd) || error("Dimensionality mismatch")
    mms = allocate_mmarrays(T, aperture_centers, maxshift)
    cms = CMStorage{T}(undef, aperture_width, maxshift; flags=flags, kwargs...)
    mismatch_apertures!(mms, fixed, moving, aperture_centers, cms, normalization=normalization)
end

"""
`mismatch_apertures!(mms, fixed, moving, aperture_centers, cms;
[normalization=:pixels])` computes the mismatch between `fixed` and
`moving` over a list of apertures at positions defined by
`aperture_centers`.  The parameters and working storage are contained
in `cms`, a `CMStorage` object. The results are stored in `mms`, an
Array-of-MismatchArrays which must have length equal to the number of
aperture centers.
"""
function mismatch_apertures!(mms, fixed, moving, aperture_centers, cms::CMStorage{T}; normalization=:pixels) where {T}
    N = ndims(cms)
    fillvalue = convert(T, NaN)
    getinds = (cms.getindices...,)::NTuple{ndims(fixed),UnitRange{Int}}
    fixedT, movingT = of_eltype(T, fixed), of_eltype(T, moving)
    for (mm, center) in zip(mms, each_point(aperture_centers))
        rng = aperture_range(center, cms.aperture_width)
        fsnip = PaddedView(fillvalue, fixedT, rng)
        erng = shiftrange.(getinds, first.(rng) .- 1)  # expanded rng
        msnip = PaddedView(fillvalue, movingT, erng)
        # Perform the calculation
        fillfixed!(cms, fsnip)
        mismatch!(mm, cms, msnip; normalization=normalization)
    end
    mms
end

# Calculate the components needed to "nancorrelate"
function fftnan!(out::NanCorrFFTs{T}, A::AbstractArray{T}, fftfunc!::Function) where {T<:Real}
    I0 = real(out.I0)
    I1 = real(out.I1)
    I2 = real(out.I2)
    _fftnan!(parent(I0), parent(I1), parent(I2), A)
    fftfunc!(out.I0)
    fftfunc!(out.I1)
    fftfunc!(out.I2)
    out
end

function _fftnan!(I0, I1, I2, A::AbstractArray{T}) where {T<:Real}
    @inbounds Threads.@threads for i in CartesianIndices(size(A))
        a = A[i]
        f = !isnan(a)
        I0[i] = f
        af = f ? a : zero(T)
        I1[i] = af
        I2[i] = af * af
    end
end

function fillfixed!(cms::CMStorage{T}, fixed::AbstractArray) where {T}
    fill!(cms.padded, NaN)
    pinds = CartesianIndices(ntuple(d -> (1:size(fixed, d)) .+ cms.maxshift[d], ndims(fixed)))
    copyto!(cms.padded, pinds, fixed, CartesianIndices(fixed))
    fftnan!(cms.fixed, cms.padded, cms.fftfunc!)
end

#### Utilities

Base.isnan(A::Array{Complex{T}}) where {T} = isnan(real(A)) | isnan(imag(A))
function sumsq_finite(A)
    s = 0.0
    for a in A
        if isfinite(a)
            s += a * a
        end
    end
    if s == 0
        error("No finite values available")
    end
    s
end

### Deprecations

function CMStorage{T}(::UndefInitializer, aperture_width::WidthLike, maxshift::WidthLike; kwargs...) where {T<:Real}
    Base.depwarn("CMStorage with aperture_width::$(typeof(aperture_width)) and maxshift::$(typeof(maxshift)) is deprecated, use tuples instead", :CMStorage)
    (N = length(aperture_width)) == length(maxshift) || error("Dimensionality mismatch")
    return CMStorage{T,N}(undef, (aperture_width...,), (maxshift...,); kwargs...)
end

@deprecate CMStorage(::Type{T}, aperture_width, maxshift; kwargs...) where {T} CMStorage{T}(undef, aperture_width, maxshift; kwargs...)

end
