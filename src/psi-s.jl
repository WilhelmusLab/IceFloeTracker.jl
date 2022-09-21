
"""
    dx, dy = grad(x::Vector{<:Number}, y::Vector{<:Number})

Make gradient vector field for the set of points with coordinates in vectors `x` and `y`. Return a tuple with `dx` and `dy` in that order. 
"""
function grad(x::Vector{<:Number}, y::Vector{<:Number})
    dx = x[2:end] - x[1:end-1]
    dy = y[2:end] - y[1:end-1]
    return dx, dy
end

"""
    dx, dy = grad(A::Matrix{<:Number})

Make gradient vector field for the set of points with coordinates in the rows of the matrix `A` with x-coordinates down column 1 and y-coordinates down column 2. Return a tuple with `dx` and `dy` in that order. 
"""
function grad(A::Matrix{<:Number})
    # Grab each col of A
    x,y = A[:,1], A[:,2]
    return grad(x,y)    
end

"""
    arclength(x::Vector{<:Number}, y::Vector{<:Number})

Return the traversed arclength along the curve represented by the points with coordinates `x` and `y`.
"""
function arclength(x::Vector{<:Number}, y::Vector{<:Number})::Vector{<:Number}
    dx,dy = grad(x,y)
    return cumsum([norm(collect(v_i)) for v_i in eachrow([dx dy])])
end

"""
    norm(v)

Get the euclidean norm of the vector `v`.
"""
function norm(v::Vector{<:Number})
    return sum(v.^2)^.5
end

"""
    atan2(y,x)

Wrapper of `Base.atan` that returns the angle of vector (x,y) in the range [0, 2π).
"""
function atan2(y::Number,x::Number)
    ang = atan(y,x)
    if y<0
        ang += 2*pi
    end
    return ang
end

"""
    make_psi_s(xs::Vector{Float64}, ys::Vector{Float64}; rangeout::Int64=0,unwrap::Bool=true, s::Bool=false)

Builds the ψ-s curve defined by vectors `xs` and `ys`.

Following the convention in [1], the ψ-s curve has values in [0, 2π) by default; use `range` to control this behavior.

See also [`bwtraceboundary`](@ref), [`resample_boundary`](@ref)

# Arguments
- `xs`: Vector of x-coordinates
- `ys`: corresponding vector of y-coordinates
- `rangeout`: 0 (default) for phase values in [0, 2π); 1 for phase values in (-π, π].
- `unwrap`: set to `true` to get "unwrapped" phases (default).
- `s`: set to `true` to additionally return the cummulative arclength traversed at each point (set `false` by default). 

# Reference
[1] McConnell, Ross, et al. "psi-s correlation and dynamic time warping: two methods for tracking ice floes in SAR images." IEEE Transactions on Geoscience and Remote sensing 29.6 (1991): 1004-1012.

# Example

The example below builds a cardioid and obtains its ψ-s curve.

```jldoctest; setup = :(using IceFloeTracker, Plots)
julia> t = range(0,2pi,201);

julia> x = @. cos(t)*(1-cos(t));

julia> y = @. sin(t)*(1-cos(t));

julia> plot(x,y) # visualize the cardiod

julia> psi = make_psi_s(x,y)
200-element Vector{Float64}:
 0.031415926535897934
 0.07330344575873542
 0.11937977657362882
 0.16605452665067502
 0.21292875044763676
 0.2598936441351966
 ⋮
 9.211849210321741
 9.258723434118705
 9.305398184195749
 9.351474515010644
 9.393362034233482

 julia> plot(t[1:end-1], psi) # inspect psi-s curve -- should be a straigth line from (0, 0) to (2π, 3π)
```
"""
function make_psi_s(x::Vector{<:Number},
                    y::Vector{<:Number};rangeout::Bool=true,
                    unwrap::Bool=true)
    # gradient
        dx, dy = grad(x, y)
        
    # get phase curve
        if rangeout
            phase = atan2.(dy, dx)
        else
            phase = atan.(dy, dx)
        end

        if unwrap
            phase = DSP.unwrap(phase)
        end
    
    # compute arclength
        s = cumsum([norm(collect(v_i)) for v_i in eachrow([dx dy])])
    
    return phase, s

end

"""
    make_psi_s(XY::Matrix{<:Number};rangeout::Bool=true,
    unwrap::Bool=true)

Alternate method of `make_psi_s` accepting input vectors `x` and `y` as a 2-column matrix `[x y]` in order to facillitate workflow (output from `resample_boundary`).
"""
function make_psi_s(XY::Matrix{<:Number};rangeout::Bool=true, unwrap::Bool=true)
    x = XY[:,1]; y=XY[:,2]
    return make_psi_s(x,y,rangeout=rangeout,unwrap=unwrap)
end
