
function save_intermediate_images(
    directory::AbstractString, results::Dict; extension::AbstractString=".png"
)
    mkdir(directory)
    for (name, image) in results
        try
            save(joinpath(directory, String(name) * extension), image)
        catch e
            @warn "an unexpected error occured saving $name"
        end
    end
    return nothing
end