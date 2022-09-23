
import DSP: xcorr

"""
    crosscorr(u::Vector{<:Number},v::Vector{<:Number};normalize::Bool=false)

Wrapper of DSP.xcorr with normalization (see https://docs.juliadsp.org/stable/convolutions/#DSP.xcorr)

"""
function crosscorr(u::Vector{<:Number},v::Vector{<:Number};normalize::Bool=false)
    c = DSP.xcorr(u,v, padmode = :longest) # same as matlab 
    if normalize
        return c / sqrt(sum(x .* x) * sum(y .* y))
    else
        return c
    end
end

