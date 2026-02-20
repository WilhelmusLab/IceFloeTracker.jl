module Preprocessing

export make_landmask_se,
    create_landmask,
    create_coastal_buffer,
    apply_landmask,
    apply_landmask!,
    LopezAcostaCloudMask,
    AbstractCloudMaskAlgorithm,
    Watkins2025CloudMask,
    create_cloudmask,
    apply_cloudmask,
    apply_cloudmask!

include("cloudmask.jl")
include("landmask.jl")

end
