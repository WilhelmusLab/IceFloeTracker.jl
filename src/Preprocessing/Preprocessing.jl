module Preprocessing

export make_landmask_se,
    create_landmask,
    create_coastal_buffer_mask,
    LopezAcostaCloudMask,
    AbstractCloudMaskAlgorithm,
    Watkins2025CloudMask,
    create_cloudmask

import ..ImageUtils: apply_mask, apply_mask!
export apply_mask, apply_mask!

include("cloudmask.jl")
include("landmask.jl")

end
