module Morphology

export bridge,
    branch,
    bwareamaxfilt,
    bwperim,
    fill_holes,
    hbreak,
    hbreak!,
    make_hbreak_dict,
    morph_fill,
    imextendedmin,
    imregionalmin,
    se_disk50,
    se_disk4,
    se_disk20,
    se_disk2,
    reconstruct,
    impose_minima,
    get_areas,
    get_max_label

include("branch.jl")
include("bridge.jl")
include("bwareamaxfilt.jl")
include("bwperim.jl")
include("fill-holes.jl")
include("hbreak.jl")
include("lut/lutbridge.jl")
include("lut/lutfill.jl")
include("minima-transform.jl")
include("morph-fill.jl")
include("operator.jl")
include("reconstruction.jl")
include("special-strels.jl")

end
