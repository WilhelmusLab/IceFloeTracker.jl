module Pipeline

export IceFloeSegmentationAlgorithm,
    LopezAcosta2019,
    LopezAcosta2019Tiling

include("abstract-algorithms.jl")
include("LopezAcosta2019.jl")
include("LopezAcosta2019Tiling.jl")

end
