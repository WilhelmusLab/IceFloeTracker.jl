# dmw: look for ways to avoid using these functions
function to_uint8(arr::AbstractMatrix{T}) where {T<:AbstractFloat}
    img = Int.(round.(arr, RoundNearestTiesAway))
    img = clamp.(img, 0, 255)
    return img
end

function to_uint8(arr::AbstractMatrix{T}) where {T<:Integer}
    img = clamp.(arr, 0, 255)
    return img
end

function to_uint8(num::T) where {T<:Union{AbstractFloat,Int,Signed}}
    num = Int(round(num, RoundNearestTiesAway))
    return clamp(num, 0, 255)
end
