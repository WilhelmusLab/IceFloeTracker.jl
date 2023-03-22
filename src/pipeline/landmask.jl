"""
    landmask(; input, output)

Given an input directory with a landmask file, create a land/soft ice mask object with both dilated and non_dilated versions. The object is serialized to the snakemake output directory. 

# Arguments
- `input`: path to image dir containing truecolor and landmask source images
- `output`: path to output dir where land-masked truecolor images and the generated binary land mask are saved
- `landmask_fname`: name of the landmask file in `input`. Default is `"landmask.tiff"`
- `outfile`: name of the serialized landmask object. Default is `"generated_landmask.jls"`
"""
function landmask(; input::String, output::String, landmask_fname::String="landmask.tiff", outfile="generated_landmask.jls")
    @info "Looking for $landmask_fname in $input"

    lmpath = joinpath(input, landmask_fname)
    check_landmask_path(lmpath)
    @info "$landmask_fname found in $input. Creating landmask..."

    img = load(lmpath)
    mkpath(output)

    # create landmask, both dilated and non-dilated as namedtuple
    serialize(joinpath(output, outfile), create_landmask(img))
    @info "Landmask created and serialized succesfully."
    return nothing
end