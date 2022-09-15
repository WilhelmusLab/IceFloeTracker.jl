
"""
    atan2(y,x)

Wrapper of `Base.atan` that returns the angle of vector (x,y) in the range [0, 2π).
"""
function atan2(y,x)
    ang = atan(y,x)
    if y<0
        return 2*pi+ang
    end
    return ang
end

function make_psi_s(xs::Vector{Float64}, ys::Vector{Float64})
    # gradient
        dx = xs[2:end] - xs[1:end-1]
        dy = ys[2:end] - ys[1:end-1]

    #angle of tangents
        angle_ = atan2.(dy,dx)
        
    # return unwrapped angle_
        return DSP.unwrap(angle_,range=2pi)
end