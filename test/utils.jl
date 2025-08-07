"""
    save_intermediate_images(
        directory::AbstractString,
        results::Dict;
        extension::AbstractString=".png",
        names::Union{AbstractArray{Symbol},Nothing}=nothing,
    )

Save images from a `results` Dict to a `directory` using the names in the dictionary as file names.
Optional arguments:
- `names` – only dictionary entries with these keys will be saved; others will be ignored
- `extension` – the file extension to use (implies the image type)
"""
function save_intermediate_images(
    directory::AbstractString,
    results::Dict;
    extension::AbstractString=".png",
    names::Union{AbstractArray{Symbol},Nothing}=nothing,
)
    mkdir(directory)
    for (name, image) in results
        # only continue if either names is undefined, or the name is in the names array
        (names === nothing || name ∈ names) || continue
        try
            save(joinpath(directory, String(name) * extension), image)
        catch e
            @warn "an unexpected error occured saving $name"
        end
    end
    return nothing
end

"""
    save_results_callback(
        path;
        extension,
        names::Union{AbstractArray{Symbol},Nothing}
    )::Function

Returns a function which saves any images 
"""
function save_results_callback(
    directory::AbstractString;
    extension::AbstractString=".png",
    names::Union{AbstractArray{Symbol},Nothing}=nothing,
)
    function callback(; kwargs...)
        mkpath(directory)
        for (name, image) in kwargs
            (names === nothing || name ∈ names) || continue
            path = joinpath(directory, String(name) * extension)
            try
                save(path, image)
            catch e
                @warn "an unexpected error occured saving $name: $e"
            end
        end
    end
    return callback
end