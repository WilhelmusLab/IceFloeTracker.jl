function check_landmask_path(lmpath::String)::Nothing
    name = basename(lmpath)
    input = dirname(lmpath)
    !isfile(lmpath) && error(
        "`$(name)` not found in $input. Please ensure a coastline image file named `$name` exists in $input.",
    )
    return nothing
end

"""
    landmask(; input, output)

Given an input directory with a landmask file and possibly truecolor images, create a land/soft ice mask. The resulting images are saved to the snakemake output directory. 

# Arguments
- `input`: path to image dir containing truecolor and landmask images
- `output`: path to output dir where land-masked truecolor images are saved

"""
function landmask(; input::String, output::String)
    landmask_fname = "landmask.tiff"
    @info "Looking for $landmask_fname in $input"

    lmpath = joinpath(input, landmask_fname)
    check_landmask_path(lmpath)
    @info "$landmask_fname found in $input. Creating landmask..."

    img = load(lmpath)
    mkpath(output)
    out = @persist create_landmask(img) joinpath(output, "generated_landmask.png")
    @info "Landmask created succefully."
    return out
end
