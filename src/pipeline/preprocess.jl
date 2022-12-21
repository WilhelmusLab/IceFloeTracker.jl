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
- `input`: path to image dir containing truecolor and landmask source images
- `output`: path to output dir where land-masked truecolor images and the generated binary land mask are saved

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
    @info "Landmask created succesfully."
    return out
end

function landmask_truecolor(lm, input::String, output::String)::Nothing
    # find truecolor imgs in input dir
    tc = sort([img for img in readdir(input) if contains(img, "truecolor")])
    total_tc = length(tc)
    @info "Found $(total_tc) truecolor images in $input. Landmasking true color images..."
    @simd for img in tc
        fname = joinpath(output, img*"-landmasked.png")
        @persist apply_landmask(load(img), lm) fname
    end
    return nothing
end