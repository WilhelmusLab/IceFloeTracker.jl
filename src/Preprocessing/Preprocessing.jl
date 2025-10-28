module Preprocessing

export make_landmask_se,
    create_landmask,
    apply_landmask,
    apply_landmask!,
    LopezAcostaCloudMask,
    AbstractCloudMaskAlgorithm,
    create_cloudmask,
    apply_cloudmask,
    apply_cloudmask!,
    create_clouds_channel

include("cloudmask.jl")
include("landmask.jl")

end
