module IceFloeTracker
include(joinpath(@__DIR__, "Landmask.jl"))

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

function landmask(;metadata::AbstractString, input::AbstractString, output::AbstractString)
  landmask_image = TiffImages.load(joinpath(input,"landmask.tiff"), mmap=true)
  landmask_binary = create_landmask(landmask_image, 50, 15)
  imagepaths = filter(endswith(r"(tiff)"), readdir(joinpath(input, "truecolor");join=true))
  images = map(x->TiffImages.load(x), imagepaths)
  masked_images = map(y->apply_landmask(y, landmask_binary), images)
  mkpath("$output")
  for (filepath, img) in collect(zip(imagepaths, masked_images))
    filename = basename(filepath)
    TiffImages.save("$output/masked_$filename", img)
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
