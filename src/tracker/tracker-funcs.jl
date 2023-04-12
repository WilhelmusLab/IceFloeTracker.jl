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
Container for matched pairs of floes. `day1props` and `day2props` are dataframes with the same column names as the input dataframes. `ratios` is a dataframe with column names `area`, `majoraxis`, `minoraxis`, `convexarea`, `area_under`, and `corr` for similarity ratios.
"""
struct MatchedPairs
    day1props::DataFrame
    day2props::DataFrame
    ratios::DataFrame
    dist::Vector{Float64}
end

"""
    MatchedPairs(df)

Return an object of type `MatchedPairs` with an empty dataframe with the same column names as `df`, an empty dataframe with column names `area`, `majoraxis`, `minoraxis`, `convexarea`, `area_under`, and `corr` for similarity ratios, and an empty vector for distances.
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
    append!(match_total.day1props, matched_pairs.day1props)
    append!(match_total.day2props, matched_pairs.day2props)
    append!(match_total.ratios, matched_pairs.ratios)
    append!(match_total.dist, matched_pairs.dist)
    return nothing
end

"""
    addmatch!(matched_pairs, newmatch)

Add `newmatch` to `matched_pairs`.
"""
function addmatch!(matched_pairs::MatchedPairs, newmatch)
    push!(matched_pairs.day1props, newmatch.day1props)
    push!(matched_pairs.day2props, newmatch.day2props)
    push!(matched_pairs.ratios, newmatch.ratios)
    push!(matched_pairs.dist, newmatch.dist)
    return nothing
end

function isempty(matched_pairs::MatchedPairs)
    return isempty(matched_pairs.day1props) &&
           isempty(matched_pairs.day2props) &&
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
struct Tracked
    data::Vector{MatchedPairs}
end

"""
    sort!(tracked::Tracked)

Sort the floes in `tracked` by area in descending order.
"""
function sort!(tracked::Tracked)
    for container in tracked.data
        p = sortperm(container.day1props, "area"; rev=true)
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
    return props_day[r, [:x, :y]]
end

"""
absdiffmeanratio(x,y)

Calculate the absolute difference between `x` and `y` divided by the mean of `x` and `y`.
"""
function absdiffmeanratio(x::Float64, y::Float64)::Float64
    return abs(x - y) / mean(x, y)
end

"""
mean(x,y)

Compute the mean of `x` and `y`.
"""
function mean(x::T, y::T)::T where {T<:Float64}
    return (x + y) / 2
end

"""
    propmatrix2df(propmatrix)

Convert floe properties matrix `propmatrix` to a dataframe-like struct with fields `area`, `majoraxis`, `minoraxis`, `convexarea`, `x_coord`, `y_coord` for the tacker step. Used for development purposes. TODO: remove this function when development is complete.
"""
function propmatrix2df(propmatrix)
    # convert props_day1 to dataframe with column names area, majoraxis, minoraxis, convexarea, coords
    return DataFrame(
        propmatrix,
        [
            :area,
            :perimeter,
            :majoraxis,
            :minoraxis,
            :orientation,
            :x,
            :y,
            :convexarea,
            :solidity,
            :left,
            :top,
            :width,
            :height,
        ],
    )
    # Props((col for col in eachcol(propmatrix[:, [1, 3, 4, 8, 6, 7]]))...)
end

"""
    poprow!(df::DataFrame, idx::Int64)

"""
function poprow!(df::DataFrame, idx::Int64)
    deleteat!(df, idx)
    return nothing
end

"""
    makeemptydffrom(df::DataFrame)

Return an object with an empty dataframe with the same column names as `df` and an empty dataframe with column names `area`, `majoraxis`, `minoraxis`, `convexarea`, `area_under`, and `corr` for similarity ratios. 
"""
function makeemptydffrom(df::DataFrame)
    return MatchingProps(
        Vector{Int64}(), similar(df, 0), makeemptyratiosdf(), Vector{Float64}()
    )
end

"""
    makeemptyratiosdf()

Return an empty dataframe with column names `area`, `majoraxis`, `minoraxis`, `convexarea`, `area_under`, and `corr` for similarity ratios.
"""
function makeemptyratiosdf()
    return DataFrame(;
        area=Float64[],
        majoraxis=Float64[],
        minoraxis=Float64[],
        convexarea=Float64[],
        area_under=Float64[],
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
    trackercond2(area1, ratios, t2=(area=1200, arearatio=0.28, majaxisratio=0.10, minaxisratio=0.12, convexarearatio=0.14))

Return `true` if the floe at `p1` and the floe at `p2` are within a certain distance of each other and the displacement time is within a certain range. Return `false` otherwise.
"""
function trackercond2(
    area1,
    ratios,
    t2=(
        area=1200,
        arearatio=0.28,
        majaxisratio=0.10,
        minaxisratio=0.12,
        convexarearatio=0.14,
    ),
)
    return area1 > t2.area &&
           ratios.area < t2.arearatio &&
           ratios.majoraxis < t2.majaxisratio &&
           ratios.minoraxis < t2.minaxisratio &&
           ratios.convexarea < t2.convexarearatio
end

"""
    trackercond3(area1, ratios, t3=(area=1200, arearatio=0.18, majaxisratio=0.07, minaxisratio=0.08, convexarearatio=0.09))

Return `true` if the floe at `p1` and the floe at `p2` are within a certain distance of each other and the displacement time is within a certain range. Return `false` otherwise.
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
    return area1 < t3.area &&
           ratios.area < t3.arearatio &&
           ratios.majoraxis < t3.majaxisratio &&
           ratios.minoraxis < t3.minaxisratio &&
           ratios.convexarea < t3.convexarearatio
end

"""
    callmatchcorr(conditions)

Condition to decide whether match_corr should be called.
"""
function callmatchcorr(conditions)
    return conditions.cond1 && (conditions.cond2 || conditions.cond3)
end

"""
    isfloegoodmatch(conditions, mct, area_under, corr)

Return `true` if the floes are a good match as per the set thresholds. Return `false` otherwise.

# Arguments
- `conditions`: tuple of booleans for evaluating the conditions
- `mct`: tuple of thresholds for the match correlation test
- `area_under` and `corr`: values returned by `match_corr`
"""
function isfloegoodmatch(conditions, mct, area_under, corr)
    return (
        (conditions.cond3 && area_under < mct.area3) ||
        (conditions.cond2 && area_under < mct.area2)
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
    majoraxisratio = absdiffmeanratio(props_day1.majoraxis[r], props_day2.majoraxis[s])
    minoraxisratio = absdiffmeanratio(props_day1.minoraxis[r], props_day2.minoraxis[s])
    convexarearatio = absdiffmeanratio(props_day1.convexarea[r], props_day2.convexarea[s])
    return (
        area=arearatio,
        majoraxis=majoraxisratio,
        minoraxis=minoraxisratio,
        convexarea=convexarearatio,
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
    return sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)
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
    return properties[dayidx], properties[dayidx + 1]
end

"""
    addmorematches!(matchingfloes, props_day, numtodadd=5)

Add `numtodadd` more matches to `matchingfloes` using the floes in `props_day`.
"""
function addmorematches!(matchingfloes, props_day, numtodadd=5)
    for j in 1:numtodadd
        s = rand(1:4)
        appendrow!(matchingfloes, props_day[s, :], randratios(), s)
    end
    return nothing
end

"""
    getbestmatchdata(idx, r, props_day1, matching_floes)

Collect the data for the best match between the `r`th floe in `props_day1` and the `idx`th floe in `matching_floes`. Return a tuple of the floe properties for day 1 and day 2 and the ratios.
"""
function getbestmatchdata(idx, r, props_day1, matching_floes)
    return (
        day1props=props_day1[r, :],
        day2props=matching_floes.props[idx, :],
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
    return filter(r -> r.x1 != 0, collisions)[:, 1:(end - 1)]
end

function deletematched!(
    (propsday1, propsday2)::Tuple{DataFrame,DataFrame}, matched::MatchedPairs
)
    deletematched!(propsday1, matched.day1props)
    deletematched!(propsday2, matched.day2props)
    return nothing
end

"""
    deleteallbut!(matched_pairs, idxs, keeper)

Delete all rows in `matched_pairs` except for the row with index `keeper` in `idxs`.
"""
function deleteallbut!(matched_pairs, idxs, keeper)
    for i in idxs
        if i !== keeper
            deleteat!(matched_pairs.ratios, i)
            deleteat!(matched_pairs.day1props, i)
            deleteat!(matched_pairs.day2props, i)
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
    collisions = getcollisionslocs(matched_pairs.day2props)
    for collision in collisions
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

# generate random ratio data
function randratios()
    return (
        area=rand(),
        majoraxis=rand(),
        minoraxis=rand(),
        convexarea=rand(),
        area_under=rand(),
        corr=rand(),
    )
end

# match_corr related functions
"""
    TODO
    matchcorr(r, s, dayi, delta_time, floe_library)

Compute the match correlation between the `r`th floe in `props_day1` and the `s`th floe in `props_day2`. Return a tuple of the area under the curve and their correlation score.
"""
function matchcorr(r, s, dayi, delta_time, floe_library)
    return (area_under=rand() / 4, corr=1 - rand() / 3)
end

"""
    corr(f1,f2)

Return the correlation between the psi-s curves `p1` and `p2`.
"""
function corr(p1, p2)
    cc, _ = maximum.(IceFloeTracker.crosscorr(p1, p2; normalize=true))
    return cc
end

"""
    myargmax(mms)

Return the index of the maximum value in `mms` and the maximum value.
"""
myargmax(mms) = (max=maximum(mms), idx=argmax(mms) + 1)

myargmin(mms) = (min=minimum(mms), idx=argmin(mms) + 1)

"""
    computethetarevROT180(revised,t=180)

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

function addψs!(props)
    props.psi = map(buildψs, props.mask)
    return nothing
end
