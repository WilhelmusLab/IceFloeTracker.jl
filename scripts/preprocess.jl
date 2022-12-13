using IceFloeTracker: create_landmask, make_landmask_se, load, Gray, @persist, writedlm
imshow(x)=Gray.(x)
landmask_file = "landmask.tiff"
input = joinpath(@__DIR__,"input")
output= joinpath(@__DIR__,"$output")

# Check input and output dirs. TODO: refacor as function
!isdir(input) && throw(Base.IOError("$input", 999))
!isdir(output) && throw(Base.IOError("$output", 999))
files = readdir(input)

# Check dist of files TODO: refactor as function
aquaimgs = [f for f in files if contains(f,"aqua")]
terraimgs = [f for f in files if contains(f,"terra")]
countaqua != countterra && error("Count of aqua/terra images do not match:\ncount aqua = $countaqua \ncount terra = $countterra")

truecolor = [f for f in files if contains(f,"truecolor")]; sort!(truecolor)
reflectance = [f for f in files if contains(f,"reflectance")]; sort!(reflectance)

length(truecolor) != length(reflectance) && error("Count of truecolor/reflectance images do not match:\ncount truecolor = $counttruecolor \ncount reflectance = $countreflectance")

tc_ref_pairs = zip(truecolor, reflectance)

function check_landmask_path(lmpath::String)::Nothing
    !isfile(lmpath) && error("`$(landmask_file)` not found in $dir. Please ensure a 
    coastline image file named `$(landmask_file)` exists in $input.")
end

function landmask(; metadata::String,
    landmask_fname::String="landmask.tiff",
    input_dir::String, output_dir::String)

    @info "Looking for $landmask_fname in $input_dir"
    lmpath = joinpath(input_dir,landmask_fname)
    check_landmask_path(lmpath)
    

    @info "$landmask_file found in $input_dir. Creating landmask..."

    img = load(lmpath)
    landmask = create_landmask(img)

    @persist landmask joinpath(output_dir, "landmask.png")

    @info "Landmask created succefully."
    return landmask
end
