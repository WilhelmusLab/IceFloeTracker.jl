module IceFloeTracker

# Includes in order of dependencies:
# No dependencies first, then modules with dependencies on those, etc.
# Each module should list all of its global exports at the top of its file.
include("skimage/skimage.jl")
include("Utils/Utils.jl")
include("Data/Data.jl")
include("ImageUtils/ImageUtils.jl")
include("Geospatial/Geospatial.jl")
include("Morphology/Morphology.jl")
include("Filtering/Filtering.jl")
include("Preprocessing/Preprocessing.jl")
include("Segmentation/Segmentation.jl")
include("Tracking/Tracking.jl")
include("Pipeline/LopezAcosta2019Tiling.jl")
include("Pipeline/LopezAcosta2019.jl")

end
