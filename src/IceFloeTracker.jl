module IceFloeTracker
using Images
using ImageProjectiveGeometry
using DelimitedFiles
using Dates
using ImageContrastAdjustment
using ImageSegmentation
using Peaks
using StatsBase
using Interpolations
using DataFrames
using PyCall
using Clustering
using DSP

include("utils.jl")
include("persist.jl")
include("landmask.jl")
include("cloudmask.jl")
include("normalization.jl")
include("ice-water-discrimination.jl")
include("anisotropic_image_diffusion.jl")
include("bwtraceboundary.jl")
include("resample-boundary.jl")
include("psi-s.jl")

const sk_measure = PyNULL()

function __init__()
    return copy!(sk_measure, pyimport_conda("skimage.measure", "scikit-image"))
end

include("regionprops.jl")
include("segmentation_a_direct.jl")
include("segmentation_b.jl")
include("segmentation_c.jl")
include("segmentation_d_e.jl")

function fetchdata(; output::AbstractString)
    mkpath("$output")
    touch("$output/metadata.json")

    mkpath("$output/landmask")
    touch("$output/landmask/landmask.tiff")

    mkpath("$output/truecolor")
    touch("$output/truecolor/a.tiff")
    touch("$output/truecolor/b.tiff")
    touch("$output/truecolor/c.tiff")

    mkpath("$output/reflectance")
    touch("$output/reflectance/a.tiff")
    touch("$output/reflectance/b.tiff")
    touch("$output/reflectance/c.tiff")
    return nothing
end

"""
    landmask(;metadata, input, output)

Given an input directory with a landmask file and truecolor images, create a land/soft ice mask, and apply to the truecolor images. The resulting land-masked images are saved to the snakemake output directory. 

# Arguments
- `metadata`: JSON file with metadata
- `input`: path to image dir containing truecolor and landmask images
- `output`: path to output dir where land-masked truecolor images are saved

"""
function landmask(; metadata::AbstractString, input::AbstractString, output::AbstractString)
    landmask_image = TiffImages.load(joinpath(input, "landmask.tiff"); mmap=true)
    landmask_binary = create_landmask(landmask_image)
    imagepaths = filter(endswith("tiff"), readdir(joinpath(input, "truecolor"); join=true))
    mkpath("$output")

    for imagepath in imagepaths
        image = TiffImages.load(imagepath)
        image = apply_landmask(image, landmask_binary)
        filename = basename(imagepath)
        TiffImages.save("$output/masked_$filename", image)
    end
    return nothing
end

function cloudmask(;
    metadata::AbstractString, input::AbstractString, output::AbstractString
)
    mkpath("$output")
    touch("$output/a.tiff")
    touch("$output/b.tiff")
    touch("$output/c.tiff")
    return nothing
end

end
