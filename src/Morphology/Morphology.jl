module Morphology

export bridge, branch, bwareamaxfilt, fill_holes, hbreak, hbreak!, morph_fill

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

end
