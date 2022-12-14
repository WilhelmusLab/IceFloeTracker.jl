using IceFloeTracker: create_landmask, load, Gray, @persist, writedlm

function check_landmask_path(lmpath::String)::Nothing
    name = basename(lmpath)
    input = dirname(lmpath)
    !isfile(lmpath) && error("`$(name)` not found in $input. Please ensure a coastline image file named `$name` exists in $input.")
    nothing
end

function landmask(; metadata::Union{Nothing, String}=nothing,
    landmask_fname::Union{Nothing, String}=nothing,
    input_dir::String, output_dir::String)

    # use default name for landmask if not provided
    isnothing(landmask_fname) && (landmask_fname = "landmask.tiff")
    @info "Looking for $landmask_fname in $input_dir"
    
    lmpath = joinpath(input_dir,landmask_fname)
    check_landmask_path(lmpath)
    @info "$landmask_fname found in $input_dir. Creating landmask..."

    img = load(lmpath)
    mkpath(output_dir)
    @persist create_landmask(img) joinpath(output_dir, "generated_landmask.png")
    @info "Landmask created succefully."
    return nothing
end
