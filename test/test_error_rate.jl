# TODO: make this into a file not named test...
function test_similarity(imgA::BitMatrix, imgB::BitMatrix, error_rate=0.005)
    error = sum(imgA .!== imgB) / prod(size(imgA))
    res = error_rate > error
    if res
        @info "Test passed with $error mismatch with threshold $error_rate"
    else
        @warn "Test failed with $error mismatch with threshold $error_rate"
    end
    return res
end
