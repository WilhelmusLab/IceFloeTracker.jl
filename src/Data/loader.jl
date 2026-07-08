export AbstractLoader, GitHubLoader

abstract type AbstractLoader end

import Downloads: RequestError, download

@kwdef struct GitHubLoader <: AbstractLoader
    url::AbstractString
    ref::AbstractString
    cache_dir::AbstractString
end

function (p::GitHubLoader)(file::AbstractString)::Union{AbstractString,Nothing}
    source = joinpath(p.url, "raw", p.ref, file)
    target = joinpath(p.cache_dir, p.ref, file)
    data = _get_file(source, target)
    return data
end

# Number of times to attempt a download before giving up.
const _MAX_DOWNLOAD_RETRIES = 5

function _get_file(
    file_url::AbstractString,
    file_path::AbstractString;
    retries::Integer=_MAX_DOWNLOAD_RETRIES,
)::Union{AbstractString,Nothing}
    @debug "looking for file at $(file_path). File exists: $(isfile(file_path))"
    isfile(file_path) && return file_path

    mkpath(dirname(file_path))
    for attempt in 1:retries
        try
            download(file_url, file_path)
            return file_path
        catch e
            isa(e, RequestError) || rethrow(e)

            status = e.response === nothing ? 0 : e.response.status
            # A genuine "not found" is not retryable; report the file as missing.
            if status == 404
                @debug "nothing at $(file_url)"
                return nothing
            end

            # Transient failures (e.g. HTTP 429 rate limiting, 5xx server
            # errors, or network errors) are retried with exponential backoff
            # and jitter to avoid hammering the server.
            if attempt == retries
                @warn "failed to download $(file_url) after $(retries) attempts" status
                rethrow(e)
            end
            backoff = 2.0^(attempt - 1) + rand()
            @debug "retrying download of $(file_url) (attempt $(attempt)/$(retries), status $(status)); sleeping $(round(backoff; digits=2))s"
            sleep(backoff)
        end
    end
end
