
export AbstractLoader, GitHubLoader

abstract type AbstractLoader end

import Downloads: RequestError, download

@kwdef struct GitHubLoader <: AbstractLoader
    url::AbstractString
    ref::AbstractString
    cache_dir::AbstractString
end

function (p::GitHubLoader)(file::AbstractString)::AbstractString
    source = joinpath(p.url, "raw", p.ref, file)
    target = joinpath(p.cache_dir, p.ref, file)
    data = _get_file(source, target)
    return data
end

function _get_file(file_url::AbstractString, file_path::AbstractString)::AbstractString
    @debug "looking for file at $(file_path). File exists: $(isfile(file_path))"
    if !isfile(file_path)
        try
            mkpath(dirname(file_path))
            download(file_url, file_path)
        catch e
            if isa(e, RequestError)
                @debug "nothing at $(file_url)"
                return nothing
            else
                rethrow(e)
            end
        end
    end
    return file_path
end