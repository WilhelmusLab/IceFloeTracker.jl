module IceFloeTracker
include(joinpath(@__DIR__, "landmask.jl"))

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
    landmask(metadata, input, output)

Call the `create_landmask` and `apply_landmask` functions to zero out pixels over land in truecolor images, including the soft ice region adjacent to land, and return land-masked images.

# Arguments
- `metadata`: JSON file with metadata
- `input`: path to image dir containing truecolor and landmask images
- `output`: path to function output dir

"""
function landmask(;metadata::AbstractString, input::AbstractString, output::AbstractString)
  landmask_image = TiffImages.load(joinpath(input,"landmask.tiff"), mmap=true)
  landmask_binary = create_landmask(landmask_image, 50, 15)
  imagepaths = filter(endswith("tiff"), readdir(joinpath(input, "truecolor");join=true))
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
