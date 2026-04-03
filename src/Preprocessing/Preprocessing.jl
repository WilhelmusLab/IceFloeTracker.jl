module Preprocessing

export make_landmask_se,
    create_landmask,
    LopezAcostaCloudMask,
    AbstractCloudMaskAlgorithm,
    Watkins2025CloudMask,
    create_cloudmask

import ..ImageUtils: apply_mask, apply_mask!
export apply_mask, apply_mask!

include("cloudmask.jl")
include("landmask.jl")

end
