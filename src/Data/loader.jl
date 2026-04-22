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

function _get_file(
    file_url::AbstractString,
    file_path::AbstractString;
    max_attempts::Integer=3,
    download_fn::Function=download,
    validate_fn::Function=load,
)::AbstractString
    mkpath(dirname(file_path))

    if isfile(file_path)
        try
            validate_fn(file_path)
            return file_path
        catch e
            @info "file validation failed, removing $(file_path)" exception = e
            rm(file_path; force=true)
        end
    end

    for attempt in 1:max_attempts
        @show attempt
        try
            download_fn(file_url, file_path)
            validate_fn(file_path)
            return file_path
        catch e
            if isa(e, RequestError) && e.code != 429
                @debug "download failed for $(file_url) on attempt $(attempt)" exception = e
                continue
            else
                rethrow(e)
            end
        end
    end
    throw(
        ErrorException("Failed to get file from $(file_url) after $(max_attempts) attempts")
    )
end
