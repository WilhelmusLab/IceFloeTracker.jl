function imcomplement(img::Matrix{T}) where {T<:Union{Unsigned,Int}}
    return 255 .- img
end

function imcomplement(img::Matrix{Gray{Float64}})
    return 1 .- img
end
