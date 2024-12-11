# need this for adding methods to Base functions
import Base.isempty
import Base.isequal

# Containers and methods for preliminary matches
struct MatchingProps
    s::Vector{Int64}
    props::DataFrame
    ratios::DataFrame
    dist::Vector{Float64}
end

"""
Container for matched pairs of floes. `props1` and `props2` are dataframes with the same column names as the input dataframes. `ratios` is a dataframe with column names `area`, `majoraxis`, `minoraxis`, `convex_area`, `area_mismatch`, and `corr` for similarity ratios. `dist` is a vector of (pixel) distances between paired floes.
"""
struct MatchedPairs
    props1::DataFrame
    props2::DataFrame
    ratios::DataFrame
    dist::Vector{Float64}
end

"""
    MatchedPairs(df)

Return an object of type `MatchedPairs` with an empty dataframe with the same column names as `df`, an empty dataframe with column names `area`, `majoraxis`, `minoraxis`, `convex_area`, `area_mismatch`, and `corr` for similarity ratios, and an empty vector for distances.
"""
function MatchedPairs(df)
    emptypropsdf = similar(df, 0)
    return MatchedPairs(emptypropsdf, copy(emptypropsdf), makeemptyratiosdf(), Float64[])
end

"""
    update!(match_total::MatchedPairs, matched_pairs::MatchedPairs)

Update `match_total` with the data from `matched_pairs`.
"""
function update!(match_total::MatchedPairs, matched_pairs::MatchedPairs)
    append!(match_total.props1, matched_pairs.props1)
    append!(match_total.props2, matched_pairs.props2)
    append!(match_total.ratios, matched_pairs.ratios)
    append!(match_total.dist, matched_pairs.dist)
    return nothing
end

"""
    addmatch!(matched_pairs, newmatch)

Add `newmatch` to `matched_pairs`.
"""
function addmatch!(matched_pairs::MatchedPairs, newmatch)
    push!(matched_pairs.props1, newmatch.props1)
    push!(matched_pairs.props2, newmatch.props2)
    push!(matched_pairs.ratios, newmatch.ratios)
    push!(matched_pairs.dist, newmatch.dist)
    return nothing
end

function isempty(matched_pairs::MatchedPairs)
    return isempty(matched_pairs.props1) &&
           isempty(matched_pairs.props2) &&
           isempty(matched_pairs.ratios)
end

"""
    appendrows!(df::MatchingProps, props::T, ratios, idx::Int64, dist::Float64) where {T<:DataFrameRow}

Append a row to `df.props` and `df.ratios` with the values of `props` and `ratios` respectively.
"""
function appendrows!(
    df::MatchingProps, props::T, ratios, idx::Int64, dist::Float64
) where {T<:DataFrameRow}
    push!(df.s, idx)
    push!(df.props, props)
    push!(df.ratios, ratios)
    push!(df.dist, dist)
    return nothing
end

"""
    isequal(matchedpairs1::MatchedPairs, matchedpairs2::MatchedPairs)

Return `true` if `matchedpairs1` and `matchedpairs2` are equal, `false` otherwise.
"""
function isequal(matchedpairs1::MatchedPairs, matchedpairs2::MatchedPairs)
    return all((
        isequal(getfield(matchedpairs1, name), getfield(matchedpairs2, name)) for
        name in namesof(matchedpairs1)
    ))
end

# Final pairs container and associated methods
"""

Container for final matched pairs of floes. `data` is a vector of `MatchedPairs` objects.
"""
struct Tracked
    data::Vector{MatchedPairs}
end

"""
    sort!(tracked::Tracked)

Sort the floes in `tracked` by area in descending order.
"""
function Base.sort!(tracked::Tracked)
    for container in tracked.data
        p = sortperm(container.props1, "area"; rev=true)
        for nm in namesof(container)
            getproperty(container, nm)[:, :] = getproperty(container, nm)[p, :]
        end
    end
    return nothing
end

function Tracked()
    return Tracked(MatchedPairs[])
end

function update!(tracked::Tracked, matched_pairs::MatchedPairs)
    push!(tracked.data, matched_pairs)
    return nothing
end

# Misc functions
"""
modenan(x::AbstractVector{<:Float64})

Return the mode of `x` or `NaN` if `x` is empty.
"""
function modenan(x::AbstractVector{<:Int64})
    length(x) == 0 && return NaN
    return mode(x)
end

"""
getcentroid(props_day::DataFrame, r)

Get the coordinates of the `r`th floe in `props_day`.
"""
function getcentroid(props_day::DataFrame, r)
    return props_day[r, [:row_centroid, :col_centroid]]
end

"""
absdiffmeanratio(x, y)

Calculate the absolute difference between `x` and `y` divided by the mean of `x` and `y`.
"""
function absdiffmeanratio(x::T, y::T)::Float64 where {T<:Real}
    return abs(x - y) / mean(x, y)
end

"""
mean(x,y)

Compute the mean of `x` and `y`.
"""
function mean(x::T, y::T)::Float64 where {T<:Real}
    return (x + y) / 2
end

"""
    makeemptydffrom(df::DataFrame)

Return an object with an empty dataframe with the same column names as `df` and an empty dataframe with column names `area`, `majoraxis`, `minoraxis`, `convex_area`, `area_mismatch`, and `corr` for similarity ratios.
"""
function makeemptydffrom(df::DataFrame)
    return MatchingProps(
        Vector{Int64}(), similar(df, 0), makeemptyratiosdf(), Vector{Float64}()
    )
end

"""
    makeemptyratiosdf()

Return an empty dataframe with column names `area`, `majoraxis`, `minoraxis`, `convex_area`, `area_mismatch`, and `corr` for similarity ratios.
"""
function makeemptyratiosdf()
    return DataFrame(;
        area=Float64[],
        majoraxis=Float64[],
        minoraxis=Float64[],
        convex_area=Float64[],
        area_mismatch=Float64[],
        corr=Float64[],
    )
end

#= Conditions for a match:
    Condition 1: displacement time delta =#

"""
    trackercond1(p1, p2, delta_time, t1=(dt = (30, 100, 1300), dist=(15, 30, 120)))

Return `true` if the floe at `p1` and the floe at `p2` are within a certain distance of each other and the displacement time is within a certain range. Return `false` otherwise.

# Arguments
- `p1`: coordinates of floe 1
- `p2`: coordinates of floe 2
- `delta_time`: time elapsed from image day 1 to image day 2
- `t`: tuple of thresholds for elapsed time and distance

"""
function trackercond1(d, delta_time, t1=(dt=(30, 100, 1300), dist=(15, 30, 120)))
    return (delta_time < t1.dt[1] && d < t1.dist[1]) ||
           (delta_time >= t1.dt[1] && delta_time <= t1.dt[2] && d < t1.dist[2]) ||
           (delta_time >= t1.dt[3] && d < t1.dist[3])
end

"""
    trackercond2(area1, ratios, t2=(area=1200, arearatio=0.28, majaxisratio=0.10, minaxisratio=0.12, convex_area=0.14))

Set of conditions for "big" floes. Return `true` if the area of the floe is greater than `t2.area` and the similarity ratios are less than the corresponding thresholds in `t2`. Return `false` otherwise.
"""
function trackercond2(
    area1,
    ratios,
    t2=(area=1200, arearatio=0.28, majaxisratio=0.10, minaxisratio=0.12, convex_area=0.14),
)
    return area1 > t2.area &&
           ratios.area < t2.arearatio &&
           ratios.majoraxis < t2.majaxisratio &&
           ratios.minoraxis < t2.minaxisratio &&
           ratios.convex_area < t2.convexarearatio
end

"""
    trackercond3(area1, ratios, t3=(area=1200, arearatio=0.18, majaxisratio=0.07, minaxisratio=0.08, convex_area=0.09))

Set of conditions for "small" floes. Return `true` if the area of the floe is less than `t3.area` and the similarity ratios are less than the corresponding thresholds in `t3`. Return `false` otherwise
"""
function trackercond3(
    area1,
    ratios,
    t3=(
        area=1200,
        arearatio=0.18,
        majaxisratio=0.07,
        minaxisratio=0.08,
        convexarearatio=0.09,
    ),
)
    return area1 <= t3.area &&
           ratios.area < t3.arearatio &&
           ratios.majoraxis < t3.majaxisratio &&
           ratios.minoraxis < t3.minaxisratio &&
           ratios.convex_area < t3.convexarearatio
end

"""
    callmatchcorr(conditions)

Condition to decide whether match_corr should be called.
"""
function callmatchcorr(conditions)
    return conditions.cond1 && (conditions.cond2 || conditions.cond3)
end

"""
    isfloegoodmatch(conditions, mct, area_mismatch, corr)

Return `true` if the floes are a good match as per the set thresholds. Return `false` otherwise.

# Arguments
- `conditions`: tuple of booleans for evaluating the conditions
- `mct`: tuple of thresholds for the match correlation test
- `area_mismatch` and `corr`: values returned by `match_corr`
"""
function isfloegoodmatch(conditions, mct, area_mismatch, corr)
    return (
        (conditions.cond3 && area_mismatch < mct.area3) ||
        (conditions.cond2 && area_mismatch < mct.area2)
    ) && corr > mct.corr
end

"""
    compute_ratios((props_day1, r), (props_day2,s))

Compute the ratios of the floe properties between the `r`th floe in `props_day1` and the `s`th floe in `props_day2`. Return a tuple of the ratios.

# Arguments
- `props_day1`: floe properties for day 1
- `r`: index of floe in `props_day1`
- `props_day2`: floe properties for day 2
- `s`: index of floe in `props_day2`
"""
function compute_ratios((props_day1, r), (props_day2, s))
    arearatio = absdiffmeanratio(props_day1.area[r], props_day2.area[s])
    majoraxisratio = absdiffmeanratio(
        props_day1.major_axis_length[r], props_day2.major_axis_length[s]
    )
    minoraxisratio = absdiffmeanratio(
        props_day1.minor_axis_length[r], props_day2.minor_axis_length[s]
    )
    convex_area = absdiffmeanratio(props_day1.convex_area[r], props_day2.convex_area[s])
    return (
        area=arearatio,
        majoraxis=majoraxisratio,
        minoraxis=minoraxisratio,
        convex_area=convex_area,
    )
end

"""
    compute_ratios_conditions((props_day1, r), (props_day2, s), delta_time, t)

Compute the conditions for a match between the `r`th floe in `props_day1` and the `s`th floe in `props_day2`. Return a tuple of the conditions.

# Arguments
- `props_day1`: floe properties for day 1
- `r`: index of floe in `props_day1`
- `props_day2`: floe properties for day 2
- `s`: index of floe in `props_day2`
- `delta_time`: time elapsed from image day 1 to image day 2
- `t`: tuple of thresholds for elapsed time and distance. See `pair_floes` for details.
"""
function compute_ratios_conditions((props_day1, r), (props_day2, s), delta_time, thresh)
    t1, t2, t3 = thresh
    p1 = getcentroid(props_day1, r)
    p2 = getcentroid(props_day2, s)
    d = dist(p1, p2)
    area1 = props_day1.area[r]
    ratios = compute_ratios((props_day1, r), (props_day2, s))
    cond1 = trackercond1(d, delta_time, t1)
    cond2 = trackercond2(area1, ratios, t2)
    cond3 = trackercond3(area1, ratios, t3)
    return (ratios=ratios, conditions=(cond1=cond1, cond2=cond2, cond3=cond3), dist=d)
end

"""
    dist(p1, p2)

Return the distance between the points `p1` and `p2`.
"""
function dist(p1, p2)
    return sqrt(
        (p1.row_centroid - p2.row_centroid)^2 + (p1.col_centroid - p2.col_centroid)^2
    )
end

"""
    getidxmostminimumeverything(ratiosdf)

Return the index of the row in `ratiosdf` with the most minima across its columns. If there are multiple columns with the same minimum value, return the index of the first column with the minimum value. If `ratiosdf` is empty, return `NaN`.
"""
function getidxmostminimumeverything(ratiosdf)
    nrow(ratiosdf) == 0 && return NaN
    return mode([argmin(col) for col in eachcol(ratiosdf)])
end

"""
    getpropsday1day2(properties, dayidx::Int64)

Return the floe properties for day `dayidx` and day `dayidx+1`.
"""
function getpropsday1day2(properties, dayidx::Int64)
    return copy(properties[dayidx]), copy(properties[dayidx+1])
end

"""
    getbestmatchdata(idx, r, props_day1, matching_floes)

Collect the data for the best match between the `r`th floe in `props_day1` and the `idx`th floe in `matching_floes`. Return a tuple of the floe properties for day 1 and day 2 and the ratios.
"""
function getbestmatchdata(idx, r, props_day1, matching_floes)
    matching_floes_props = matching_floes.props[idx, :]
    cols = names(matching_floes_props)
    return (
        props1=props_day1[r, cols],
        props2=matching_floes_props,
        ratios=matching_floes.ratios[idx, :],
        dist=matching_floes.dist[idx],
    )
end

"""
    getidxofrow(rw0, df)

Return the indices of the rows in `df` that are equal to `rw0`.
"""
function getidxofrow(rw0, df)
    return findall([rw0 == row for row in eachrow(df)])
end

# Collision handling

"""
    getcollisionslocs(df)

Return a vector of tuples of the row and the index of the row in `df` that has a collision with another row in `df`.
"""
function getcollisionslocs(df) #::Vector{Tuple{DataFrameRow{DataFrame, DataFrames.Index}, Vector{Int64}}}
    typeof(df) <: DataFrameRow &&
        return Tuple{DataFrameRow{DataFrame,DataFrames.Index},Vector{Int64}}[]
    collisions = getcollisions(df)
    return [(data=rw, idxs=getidxofrow(rw, df)) for rw in eachrow(collisions)]
end

namesof(obj::MatchedPairs) = fieldnames(typeof(obj))

"""
getcollisions(matchedpairs)

Get nonunique rows in `matchedpairs`.
"""
function getcollisions(matchedpairs)
    collisions = transform(matchedpairs, nonunique)
    return filter(r -> r.x1 != 0, collisions)[:, 1:(end-1)]
end

function deletematched!(
    (propsday1, propsday2)::Tuple{DataFrame,DataFrame}, matched::MatchedPairs
)
    deletematched!(propsday1, matched.props1)
    deletematched!(propsday2, matched.props2)
    return nothing
end

"""
    deleteallbut!(matched_pairs, idxs, keeper)

Delete all rows in `matched_pairs` except for the row with index `keeper` in `idxs`.
"""
function deleteallbut!(matched_pairs, idxs, keeper)
    for i in sort(idxs; rev=true)
        if i !== keeper
            deleteat!(matched_pairs.ratios, i)
            deleteat!(matched_pairs.props1, i)
            deleteat!(matched_pairs.props2, i)
            deleteat!(matched_pairs.dist, i)
        end
    end
end

"""
    ismember(df1,df2)

Return a boolean array indicating whether each row in `df1` is a member of `df2`.
"""
function ismember(df1, df2)
    return in.(eachrow(df1), Ref(eachrow(df2)))
end

function resolvecollisions!(matched_pairs)
    collisions = getcollisionslocs(matched_pairs.props2)
    for collision in reverse(collisions)
        bestentry = getidxmostminimumeverything(matched_pairs.ratios[collision.idxs, :])
        keeper = collision.idxs[bestentry]
        deleteallbut!(matched_pairs, collision.idxs, keeper)
    end
end

function deletematched!(propsday::DataFrame, matched::DataFrame)
    toremove = findall(ismember(propsday, matched))
    return deleteat!(propsday, toremove)
end

isnotnan(x) = !isnan(x)

# match_corr related functions

"""
    corr(f1,f2)

Return the correlation between the psi-s curves `p1` and `p2`.
"""
function corr(p1, p2)
    cc, _ = maximum.(IceFloeTracker.crosscorr(p1, p2; normalize=true))
    return cc
end

"""
   normalizeangle(revised,t=180)

Normalize angle to be between -180 and 180 degrees.
"""
function normalizeangle(revised, t=180)
    revised > t ? theta_revised = revised - 360 : theta_revised = revised
    return (theta_revised=theta_revised, ROT=-theta_revised)
end

function buildψs(floe)
    bd = IceFloeTracker.bwtraceboundary(floe)
    bdres = IceFloeTracker.resample_boundary(bd[1])
    return IceFloeTracker.make_psi_s(bdres)[1]
end

function addψs!(props::Vector{DataFrame})
    for prop in props
        prop.psi = map(buildψs, prop.mask)
    end
    return nothing
end

function addfloemasks!(props::Vector{DataFrame}, imgs::Vector{<:FloeLabelsImage})
    for (img, prop) in zip(imgs, props)
        IceFloeTracker.addfloemasks!(prop, img)
    end
    return nothing
end

"""
    get_unmatched(props, matched)

Return the floes in `props` that are not in `matched`.
"""
function get_unmatched(props, matched)
    _on = mapreduce(df -> Set(names(df)), intersect, [props, matched]) |> collect
    unmatched = antijoin(props, matched, on=_on)

    # Add missing columns for joining
    add_missing = ["area_mismatch", "corr"]
    [unmatched[!, n] = [missing for _ in 1:nrow(unmatched)] for n in add_missing]

    return unmatched
end

"""
    get_trajectory_heads(pairs)

Return the last row (most recent member) of each group (trajectory) in `pairs` as a dataframe.

This is used for getting the initial floe properties for the next day in search for new pairs.
"""
function get_trajectory_heads(pairs::T) where {T<:AbstractDataFrame}
    gdf = groupby(pairs, :uuid)
    return combine(gdf, last)[:, names(pairs)]
end

"""
    _swap_last_values!(df)

Swap the last two values of the `area_mismatch` and `corr` columns for each group in `df`. For bookkeeping purposes for goodness of fit data during the tracking process.
"""
function _swap_last_values!(df)
    grouped = groupby(df, :uuid)  # Group by uuid
    for sdf in grouped
        n = nrow(sdf)
        if n > 1
            # Swap last two rows for area_mismatch and corr
            sdf.area_mismatch[n], sdf.area_mismatch[n-1] = sdf.area_mismatch[n-1], sdf.area_mismatch[n]
            sdf.corr[n], sdf.corr[n-1] = sdf.corr[n-1], sdf.corr[n]
        end
    end
    return df  # The original DataFrame is modified in-place
end

"""
    get_dt(props1, r, props2, s)

Return the time difference between the `r`th floe in `props1` and the `s`th floe in `props2` in minutes.
"""
function get_dt(props1, r, props2, s)
    return (props2.passtime[s] - props1.passtime[r]) / Minute(1)
end

"""
    adduuid!(props)

Assign a unique ID to each floe in each table of floe properties.
"""
function adduuid!(props::Vector{DataFrame})
    # Assign a unique ID to each floe in each image
    for (i, prop) in enumerate(props)
        props[i].uuid = [randstring(12) for _ in 1:nrow(prop)]
    end
    return nothing
end

"""
    reset_id!(df, col)

Reset the distinct values in the column `col` of `df` to be consecutive integers starting from 1.
"""
function reset_id!(df::AbstractDataFrame, col::Union{Symbol,AbstractString}=:uuid)
    ids = unique(df[!, col])
    _map = Dict(ids .=> 1:length(ids))
    transform!(df, col => ByRow(x -> _map[x]) => col)
    return nothing
end

"""
    consolidate_matched_pairs(matched_pairs::MatchedPairs)

Consolidate the floe properties and similarity ratios of the matched pairs in `matched_pairs` into a single dataframe. Return the consolidated dataframe. Used in iteration `0`.
"""
function consolidate_matched_pairs(matched_pairs::MatchedPairs)
    # Ensure UUIDs are consistent
    matched_pairs.props2.uuid = matched_pairs.props1.uuid

    # Define columns for goodness ratios
    goodness_cols = [:area_mismatch, :corr]

    # Create top DataFrame with properties and goodness ratios
    top_df = hcat(matched_pairs.props1, matched_pairs.ratios[:, goodness_cols], makeunique=true)

    # Create missing ratios DataFrame
    missing_ratios = similar(matched_pairs.ratios[:, goodness_cols])
    missing_ratios[!, :] .= missing

    bottom_df = hcat(matched_pairs.props2, missing_ratios, makeunique=true)

    combined_df = vcat(top_df, bottom_df)

    DataFrames.sort!(combined_df, [:uuid, :passtime])

    return combined_df
end

"""
    get_matches(matched_pairs)

Return a dataframe with the properties and goodness ratios of the matched pairs (right-hand matches) in `matched_pairs`. Used in iterations `1:end`.
"""
function get_matches(matched_pairs::MatchedPairs)
    # Ensure UUIDs are consistent
    matched_pairs.props2.uuid = matched_pairs.props1.uuid

    # Define columns for goodness ratios
    goodness_cols = [:area_mismatch, :corr]

    # Create DataFrame with properties and goodness ratios
    combined_df = hcat(matched_pairs.props2, matched_pairs.ratios[:, goodness_cols], makeunique=true)

    DataFrames.sort!(combined_df, [:uuid, :passtime])

    return combined_df
end

## LatLon functions originally from IFTPipeline.jl

"""
    convertcentroid!(propdf, latlondata, colstodrop)

Convert the centroid coordinates from row and column to latitude and longitude dropping unwanted columns specified in `colstodrop` for the output data structure. Addionally, add columns `x` and `y` with the pixel coordinates of the centroid.
"""
function convertcentroid!(propdf, latlondata, colstodrop)
    latitude, longitude = [
        [latlondata[c][Int(round(x)), Int(round(y))] for
         (x, y) in zip(propdf.row_centroid, propdf.col_centroid)] for
        c in ["latitude", "longitude"]
    ]

    x, y = [
        [latlondata[c][Int(round(z))] for z in V] for
        (c, V) in zip(["Y", "X"], [propdf.row_centroid, propdf.col_centroid])
    ]

    propdf.latitude = latitude
    propdf.longitude = longitude
    propdf.x = x
    propdf.y = y
    dropcols!(propdf, colstodrop)
    return nothing
end

"""
    converttounits!(propdf, latlondata, colstodrop)

Convert the floe properties from pixels to kilometers and square kilometers where appropiate. Also drop the columns specified in `colstodrop`.
"""
function converttounits!(propdf, latlondata, colstodrop)
    if nrow(propdf) == 0
        dropcols!(propdf, colstodrop)
        insertcols!(propdf, :latitude => Float64, :longitude => Float64, :x => Float64, :y => Float64)
        return nothing
    end
    convertcentroid!(propdf, latlondata, colstodrop)
    x = latlondata["X"]
    dx = abs(x[2] - x[1])
    convertarea(area) = area * dx^2 / 1e6
    convertlength(length) = length * dx / 1e3
    propdf.area .= convertarea(propdf.area)
    propdf.convex_area .= convertarea(propdf.convex_area)
    propdf.minor_axis_length .= convertlength(propdf.minor_axis_length)
    propdf.major_axis_length .= convertlength(propdf.major_axis_length)
    propdf.perimeter .= convertlength(propdf.perimeter)
    return nothing
end

function dropcols!(df, colstodrop)
    select!(df, Not(colstodrop))
    return nothing
end
