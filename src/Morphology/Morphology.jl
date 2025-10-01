module Morphology

export bridge,
    branch,
    bwareamaxfilt,
    fill_holes,
    hbreak,
    hbreak!,
    morph_fill,
    se_disk50,
    se_disk4,
    se_disk20,
    se_disk2,
    make_landmask_se

include("lut/lutbridge.jl")
include("lut/lutfill.jl")
include("branch.jl")
include("bridge.jl")
include("bwareamaxfilt.jl")
include("bwperim.jl")
include("fill-holes.jl")
include("hbreak.jl")
include("morph-fill.jl")
include("operator.jl")
include("special-strels.jl")

end
