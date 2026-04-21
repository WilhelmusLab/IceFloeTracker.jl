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
            read(io, UInt8)
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
    max_retries::Integer=3,
    download_fn::Function=download,
    validate_fn::Function=_can_open_file,
)::AbstractString
    max_retries < 1 && throw(ArgumentError("max_retries must be at least 1."))

    @debug "looking for file at $(file_path). File exists: $(isfile(file_path))"
    mkpath(dirname(file_path))
    last_error = ErrorException("unknown download/validation failure")
    last_failure_kind = "unknown"
    for retry_attempt = 1:max_retries
        is_valid = isfile(file_path) && validate_fn(file_path)
        is_valid && return file_path

        isfile(file_path) && rm(file_path; force=true)

        try
            download_fn(file_url, file_path)
        catch e
            last_error = e
            last_failure_kind = "download"
            retry_attempt < max_retries && continue
            if isa(e, RequestError)
                @debug "failed to download $(file_url) after $(max_retries) attempts" exception = e
            end
            rethrow(e)
        end

        validate_fn(file_path) && return file_path
        last_error = ErrorException("downloaded file at $(file_path) cannot be opened")
        last_failure_kind = "validation"
    end

    throw(
        ErrorException(
            "failed to fetch valid file from $(file_url) after $(max_retries) attempts ($(last_failure_kind) failure): $(sprint(showerror, last_error))"
        ),
    )
end
