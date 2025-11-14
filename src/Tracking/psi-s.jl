import DSP: unwrap

"""
    dx, dy = grad(x::Vector{<:Number}, y::Vector{<:Number})

Make gradient vector field for the set of points with coordinates in vectors `x` and `y`. Return a tuple with `dx` and `dy` in that order. 
"""
function grad(x::Vector{<:Number}, y::Vector{<:Number})
    dx = x[2:end] - x[1:(end - 1)]
    dy = y[2:end] - y[1:(end - 1)]
    return dx, dy
end

"""
    dx, dy = grad(A::Matrix{<:Number})

Make gradient vector field for the set of points with coordinates in the rows of the matrix `A` with x-coordinates down column 1 and y-coordinates down column 2. Return a tuple with `dx` and `dy` in that order. 
"""
function grad(A::Matrix{<:Number})
    # Grab each col of A
    x, y = A[:, 1], A[:, 2]
    return grad(x, y)
end

"""
    norm(v)

Get the euclidean norm of the vector `v`.
"""
function norm(v::Vector{<:Number})
    return sum(v .^ 2)^0.5
end

"""
    atan2(y,x)

Wrapper of `Base.atan` that returns the angle of vector (x,y) in the range [0, 2π).
"""
function atan2(y::Number, x::Number)
    ang = atan(y, x)
    if y < 0
        ang += 2 * pi
    end
    return ang
end

"""
    buildψs(x::Vector{<:Number},
               y::Vector{<:Number};
               rangeout::Bool=true,
               unwrap::Bool=true)::Tuple{Vector{Float64}, Vector{Float64}}

Builds the ψ-s curve defined by vectors `x` and `y`.

Returns a tuple of vectors with the phases `ψ` in the first component and the traversed arclength in the second component. 

Following the convention in [1], the wrapped ψ-s curve has values in [0, 2π) by default; use `rangeout` to control this behavior.

See also [`bwtraceboundary`](@ref), [`resample_boundary`](@ref)

# Arguments
- `x`: Vector of x-coordinates
- `y`: corresponding vector of y-coordinates
- `rangeout`: `true` (default) for phase values in [0, 2π); `false` for phase values in (-π, π].
- `unwrap`: set to `true` to get "unwrapped" phases (default). 

# Reference
[1] McConnell, Ross, et al. "psi-s correlation and dynamic time warping: two methods for tracking ice floes in SAR images." IEEE Transactions on Geoscience and Remote sensing 29.6 (1991): 1004-1012.

# Example

The example below builds a cardioid and obtains its ψ-s curve.

```jldoctest; setup = :(using IceFloeTracker, Plots)
julia> t = range(0,2pi,201);

julia> x = @. cos(t)*(1-cos(t));

julia> y = @. sin(t)*(1-cos(t));

julia> plot(x,y) # visualize the cardioid

julia> psi, s = buildψs(x,y);

julia> [s psi] # inspect psi-s data
200×2 Matrix{Float64}:
 0.00049344  0.0314159
 0.0019736   0.0733034
 0.00444011  0.11938
 0.00789238  0.166055
 0.0123296   0.212929
 0.0177505   0.259894
 0.024154    0.306907
 0.0315383   0.35395
 0.0399017   0.401012
 0.0492421   0.448087
 ⋮
 7.96772     9.02377
 7.97511     9.07083
 7.98151     9.11787
 7.98693     9.16488
 7.99137     9.21185
 7.99482     9.25872
 7.99729     9.3054
 7.99877     9.35147
 7.99926     9.39336

 julia> plot(s, psi) # inspect psi-s curve -- should be a straight line from (0, 0) to (8, 3π)
```
"""
function buildψs(
    x::Vector{<:Number}, y::Vector{<:Number}; rangeout::Bool=true, dsp_unwrap::Bool=true
)::Tuple{Vector{Float64},Vector{Float64}}
    # gradient
    dx, dy = grad(x, y)

    # get phase curve
    if rangeout
        phase = atan2.(dy, dx)
    else
        phase = atan.(dy, dx)
    end

    if dsp_unwrap
        phase = unwrap(phase)
    end

    # compute arclength
    s = cumsum([norm(collect(v_i)) for v_i in eachrow([dx dy])])

    return phase, s
end

"""
    buildψs(XY::Matrix{<:Number};rangeout::Bool=true,
    unwrap::Bool=true)

Alternate method of `buildψs` accepting input vectors `x` and `y` as a 2-column matrix `[x y]`.
"""
function buildψs(XY::Matrix{<:Number}; rangeout::Bool=true, dsp_unwrap::Bool=true)
    x = XY[:, 1]
    y = XY[:, 2]
    return buildψs(x, y; rangeout=rangeout, dsp_unwrap=dsp_unwrap)
end


"""
    buildψs(floe_mask)

Alternate method of `buildψs` accepting binary floe mask as input.
"""
function buildψs(floe_mask::AbstractArray)
    bd = bwtraceboundary(floe_mask)
    bdres = resample_boundary(bd[1])
    return buildψs(bdres)[1]
end
