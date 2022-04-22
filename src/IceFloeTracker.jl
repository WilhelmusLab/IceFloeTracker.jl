module IceFloeTracker

function fetchdata(;output::AbstractString)
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
  mkpath("$output")
  touch("$output/a.tiff")
  touch("$output/b.tiff")
  touch("$output/c.tiff")
  return nothing
end

function cloudmask(;metadata::AbstractString, input::AbstractString, output::AbstractString)
  mkpath("$output")
  touch("$output/a.tiff")
  touch("$output/b.tiff")
  touch("$output/c.tiff")
  return nothing
end

include("landmask.jl")

end
