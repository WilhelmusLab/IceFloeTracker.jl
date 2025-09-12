function imcomplement(img::Matrix{T}) where {T<:Union{Unsigned,Int}}
    return 255 .- img
end

#TODO: Remove, and use complement.() instead
function imcomplement(img::Matrix{Gray{Float64}})
    return 1 .- img
end
