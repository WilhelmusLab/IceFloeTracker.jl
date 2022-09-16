"""
    atan2(y,x)

Wrapper of `Base.atan` that returns the angle of vector (x,y) in the range [0, 2π).
"""
function atan2(y::Vector{Float64},x::Vector{Float64})
    ang = atan.(y,x)
    for i=1:length(y)
        if y[i]<0
            ang[i] += 2*pi
        end
    end
    return ang
end

"""
    make_psi_s(xs::Vector{Float64}, ys::Vector{Float64};rangeout::Int64=0, unwrap::Bool=true)

Builds the ψ-s curve defined by vectors `xs` and `ys`.

Following the convention in [1], the ψ-s curve has values in [0, 2π) by default; use `range` to control this behavior.

See also [`bwtraceboundary`](@ref), [`resample_boundary`](@ref)

# Arguments
- `xs`: Vector of x-coordinates
- `ys`: corresponding vector of y-coordinates
- `rangeout`: 0 (default) for phase values in [0, 2π); 1 for phase values in (-π, π].
-`unwrap`: set to `true` to get "unwrapped" phases (default).

# Reference
[1] McConnell, Ross, et al. "psi-s correlation and dynamic time warping: two methods for tracking ice floes in SAR images." IEEE Transactions on Geoscience and Remote sensing 29.6 (1991): 1004-1012.

# Example

The example below builds a cardiod and obtains its ψ-s curve.

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

 julia> plot(t[1:end-1], psi) # inspect psi-s curve -- should be a straigth line from (0,0) to (2π, 3π)
```
"""
function make_psi_s(xs::Vector{Float64}, ys::Vector{Float64};rangeout::Int64=0, unwrap::Bool=true)
    @assert length(xs) == length(ys) "Vectors `xs` and `ys` must have the same size."
    @assert 0<= rangeout <= 1 "Invalid value for `rangeout`($rangeout). Choose `rangeout=0` for an phase output in [0, 2π) or `rangeout=1` for (-π, π]."

    # gradient
        dx = xs[2:end] - xs[1:end-1]
        dy = ys[2:end] - ys[1:end-1]

    # get phase curve
        if rangeout == 0
            phase = atan2(dy, dx)
        elseif rangeout == 1
            phase = atan.(dy, dx)
        end

        if unwrap
            phase = DSP.unwrap(phase)
        end
        
        return phase
end
