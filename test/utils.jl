
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