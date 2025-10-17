import DSP: xcorr
"""
    r, lags = crosscorr(u::Vector{T},
                        v::Vector{T};
                        normalize::Bool=false,
                        padmode::Symbol=:longest)

Wrapper of DSP.xcorr with normalization (see https://docs.juliadsp.org/stable/convolutions/#DSP.xcorr)

Returns the pair `(r, lags)` with the cross correlation scores `r` and corresponding `lags` according to `padmode`.

# Arguments

- `u,v`: real vectors which could have unequal length.
- `normalize`: return normalized correlation scores (`false` by default).
- `padmode`: either `:longest` (default) or `:none` to control padding of shorter vector with zeros.

# Examples

The example below builds two vectors, one a shifted version of the other, and computes various cross correlation scores.

```jldoctest; setup = :(using IceFloeTracker)
julia> n = 1:5;

julia> x = 0.48.^n;

julia> y = circshift(x,3);

julia> r, lags = crosscorr(x,y,normalize=true);

julia> [r lags]
9×2 Matrix{Float64}:
0.369648    -4.0
0.947531    -3.0
0.495695    -2.0
0.3231      -1.0
0.332519     0.0
0.15019      1.0
0.052469     2.0
0.0241435    3.0
0.00941878   4.0

julia> r, lags = crosscorr(x,y,normalize=true,padmode=:none);

julia> [r lags]
9×2 Matrix{Float64}:
0.369648    1.0
0.947531    2.0
0.495695    3.0
0.3231      4.0
0.332519    5.0
0.15019     6.0
0.052469    7.0
0.0241435   8.0
0.00941878  9.0
```

This final example builds two vectors of different length and computes some cross correlation scores.

```jldoctest; setup = :(using IceFloeTracker)
julia> n = 1:5; m = 1:3;

julia> x = 0.48.^n; y = 0.48.^m;

julia> r, lags = crosscorr(x,y,normalize=true);

julia> [r lags]
9×2 Matrix{Float64}:
0.0          -4.0
4.14728e-17  -3.0
0.178468     -2.0
0.457473     -1.0
0.994189      0.0
0.477211      1.0
0.229061      2.0
0.105402      3.0
0.0411191     4.0

julia> r, lags = crosscorr(x,y,normalize=true,padmode=:none);

julia> [r lags]
7×2 Matrix{Float64}:
0.178468   1.0
0.457473   2.0
0.994189   3.0
0.477211   4.0
0.229061   5.0
0.105402   6.0
0.0411191  7.0
```
"""
function crosscorr(
    u::Vector{T}, v::Vector{T}; normalize::Bool=false, padmode::Symbol=:longest
)::Tuple{Vector{T},Vector{T}} where {T<:Real}
    # dmw: check whether we could implement this with the tools in StatsBase or Statistics
    # so we can avoid an extra import.
    c = DSP.xcorr(u, v; padmode=padmode)
    radius = 0

    if normalize
        c = c / sqrt(sum(u .* u) * sum(v .* v))
    end

    if padmode == :longest # as in matlab 
        radius = max(length(u), length(v)) - 1
        return c, Vector((-radius):radius)
    else
        radius = length(u) + length(v) - 1
        return c, Vector(1:radius)
    end
end
