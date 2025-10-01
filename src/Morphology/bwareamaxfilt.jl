import Images: label_components

"""
    get_max_label(d::Dict{Int, Int})

Get the label `k` in dictionary `d` for which d[k] is maximal.
"""
function get_max_label(d::Dict{Int,Int})::Int
    drev = Dict{Int,Int}()
    mx = 0
    for (k, value) in d
        drev[value] = k
        if value > mx
            mx = value
        end
    end
    return drev[mx]
end

"""
    filt_except_label(labeled_arr::Array{Int64, 2}, label::Int64)

Make 0 all values in `labeled_arr` that are not equal to `label`.

See also [`filt_except_label!`](@ref) 
"""
function filt_except_label(labeled_arr::Array{Int64,2}, label::Int64)::Array{Int64,2}
    outlabel = similar(labeled_arr)
    for (i, value) in enumerate(labeled_arr)
        value != label ? outlabel[i] = 0 : outlabel[i] = label
    end
    return outlabel
end

"""
    filt_except_label!(labeled_arr::Array{Int64, 2}, label::Int64)

In-place version of `filt_except_label`.

See also [`filt_except_label`](@ref) 
"""
function filt_except_label!(labeled_arr::Array{Int64,2}, label::Int64)::Array{Int64,2}
    for (i, value) in enumerate(labeled_arr)
        value != label ? labeled_arr[i] = 0 : labeled_arr[i] = label
    end
    return labeled_arr
end

"""
    get_areas(labeled_arr::Array{T, 2})::Dict{T, Int} where T


Get the "areas" (count of pixels of a given label) of the connected components in `labeled_arr`.

Return a dictionary with the frequency distribution: label => count_of_label.
"""
function get_areas(labeled_arr::Array{T,2})::Dict{T,Int} where {T}
    d = Dict{T,Int}()
    for i in labeled_arr
        i == 0 ? continue : d[i] = get(d, i, 0) + 1
    end
    return d
end

"""
    bwareamaxfilt(bwimg::AbstractArray{Bool}, conn)

Filter the smaller (by area) connected components in `bwimg` keeping the (assumed unique) largest.

Uses 8-pixel connectivity by default (`conn=8`). Use `conn=4` for 4-pixel connectivity.

"""
function bwareamaxfilt(bwimg::AbstractArray, conn::Int=8)::BitMatrix
    return bwareamaxfilt!(copy(bwimg), conn)
end

"""
    bwareamaxfilt!(bwimg::AbstractArray)

In-place version of bwareamaxfilt.

See also [`bwareamaxfilt`](@ref) 
"""
function bwareamaxfilt!(bwimg::AbstractArray, conn::Int=8)::BitMatrix
    if conn == 8
        label = label_components(bwimg, trues(3, 3))
    elseif conn == 4
        label = label_components(bwimg)
    else
        throw(ArgumentError("Only `conn=8`(default) or `conn=4`are allowed."))
    end
    d = get_areas(label)
    mx_label = get_max_label(d)
    label = (0 .!= filt_except_label!(label, mx_label))
    for (i, value) in enumerate(bwimg)
        bwimg[i] = (value && label[i])
    end
    return bwimg
end
