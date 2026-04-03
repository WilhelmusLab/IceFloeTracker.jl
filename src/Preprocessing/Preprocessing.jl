module Preprocessing

export make_landmask_se,
    create_landmask,
    apply_mask,
    apply_mask!,
    LopezAcostaCloudMask,
    AbstractCloudMaskAlgorithm,
    Watkins2025CloudMask,
    create_cloudmask

include("cloudmask.jl")
include("landmask.jl")

end
