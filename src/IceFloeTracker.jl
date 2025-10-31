module IceFloeTracker

using Reexport

# Supporting modules
include("skimage/skimage.jl")

include("Utils/Utils.jl")
@reexport using .Utils

include("Data/Data.jl")
@reexport using .Data

include("ImageUtils/ImageUtils.jl")
@reexport using .ImageUtils

include("Geospatial/Geospatial.jl")
@reexport using .Geospatial

include("Morphology/Morphology.jl")
@reexport using .Morphology

include("Filtering/Filtering.jl")
@reexport using .Filtering

include("Preprocessing/Preprocessing.jl")
@reexport using .Preprocessing

include("Segmentation/Segmentation.jl")
@reexport using .Segmentation

include("Tracking/Tracking.jl")
@reexport using .Tracking

# Pipelines
include("Pipeline/LopezAcosta2019.jl")
export LopezAcosta2019

include("Pipeline/LopezAcosta2019Tiling.jl")
export LopezAcosta2019Tiling

end
