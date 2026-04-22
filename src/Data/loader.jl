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

function _can_open_file(file_path::AbstractString)::Bool
    try
        !isfile(file_path) && return false
        filesize(file_path) == 0 && return false
        open(file_path, "r") do io
            read(io, UInt8) # read a single UInt8 from the file
        end
        return true
    catch e
        @debug "file validation failed at $(file_path)" exception = e
        return false
    end
end

function _get_file(
    file_url::AbstractString,
    file_path::AbstractString;
    max_attempts::Integer=3,
    download_fn::Function=download,
    validate_fn::Function=_can_open_file,
)::AbstractString
    max_attempts < 1 && throw(ArgumentError("max_attempts must be at least 1."))
    mkpath(dirname(file_path))
    for attempt in 1:max_attempts
        isfile(file_path) && validate_fn(file_path) && return file_path
        isfile(file_path) && rm(file_path; force=true)
        try
            download_fn(file_url, file_path)
            validate_fn(file_path) && return file_path
        catch e
            if isa(e, RequestError) && e.code != 429
                @debug "download attempt $(attempt) failed for $(file_url)" exception = e
                attempt < max_attempts && continue
            else
                rethrow(e)
            end
        end
    end
    return error(
        "failed to fetch valid file from $(file_url) after $(max_attempts) attempts"
    )
end
