import Images: Gray

function imcomplement(img::Matrix{T}) where {T<:Union{Unsigned,Int}}
    return 255 .- img
end

# TODO: Replace all uses of this function with complement.()
function imcomplement(img::Matrix{Gray{Float64}})
    return 1 .- img
end
