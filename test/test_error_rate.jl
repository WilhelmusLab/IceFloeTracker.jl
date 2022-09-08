function test_similarity(imgA, imgB, error_rate=0.005)
    error = sum(imgA .!== imgB) / prod(size(imgA))
    res = error_rate > error
    if res
        @info "Test passed with $error mismatch with threshold $error_rate"
    else
        @warn "Test failed with $error mismatch with threshold $error_rate"
    end
    return res
end
macro test_similarity(imgA, imgB, error_rate=0.005)
    return quote
        test_similarity($(esc(imgA)), $(esc(imgB)), $(esc(error_rate)))
    end
end
