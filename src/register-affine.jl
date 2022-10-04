# Pkg.add("RegisterMismatch"); Pkg.add("RegisterQD")

function register_mismatch(fixed::AbstractArray,
                           moving::AbstractArray,
                           mxshift::Tuple(Int64,Int64)=(5,5),
                           mxrot::Float64=pi/4;
                           presmoothed::Bool=false,
                           minwidth_rot=default_minwidth_rot(fixed, SD),
                           thresh::Float64=0.1,initial_tfm=IdentityTransformation(),
                           kwargs...
                           )

    tfm, mm = qd_rigid(
                       centered(fixed), centered(moving), mxshift, mxrot,
                       presmoothed=presmoothed,
                       minwidth_rot=minwidth_rot,
                       thresh=thresh,
                       initial_tfm=initial_tfm,
                       kwargs...)
    
    return mm, tfm
end