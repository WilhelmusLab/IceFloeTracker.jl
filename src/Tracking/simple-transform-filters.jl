
export DisplacementTransform,
    TimeTransform, VelocityTransform, VelocityFilter, DisplacementFilter

function DisplacementTransform(floe)
    x0 = floe.col_centroid
    y0 = floe.row_centroid
    _inner = (x, y) -> sqrt((x - x0)^2 + (y - y0)^2)
    return [:col_centroid, :row_centroid] => ByRow(_inner) => [:Δx]
end

function TimeTransform(floe)
    _inner = (t) -> t - floe.passtime
    return :passtime => ByRow(_inner) => :Δt
end

function VelocityTransform(floe)
    _inner = (Δx, Δt) -> Δx / Δt
    return [:Δx, :Δt] => ByRow(_inner) => :u
end

function VelocityFilter(args...; umax=0.75)
    return :u => ByRow((u) -> u < (umax))
end

function DisplacementFilter(args...; umax=0.75, eps=250)
    return [:Δt, :Δx] => ByRow((Δt, Δx) -> Δx < (Δt * umax + eps))
end

function AreaFilter(args...; minimum_area=100, maximum_area=90e3)
    return :area => ByRow((a) -> a >= minimum_area && a <= maximum_area)
end

function simple_floe_tracker(
    props::Vector{DataFrame};
    candidate_transform_functions::Vector{
        Pair{Union{Symbol,Vector{Symbol}},Pair{Function,Union{Symbol,Vector{Symbol}}}}
    },
    filter_function::Vector{Pair{Union{Symbol,Vector{Symbol}},Function}},
    matching_function::Base.Callable,
    minimum_area=100,
    maximum_area=90e3,
    maximum_time_step=Day(2),
)

    # dmw: give users option to copy props rather than modify in place?
    # props = subset.(props, AreaFilter(; minimum_area, maximum_area))

    # Start_new_trajectory adds head_uuid and trajectory_uuid columns to props
    # The starting trajectories are just the floes visible and large enough on day 1.
    init_idx = 1
    nrow(props[init_idx]) == 0 && begin
        while nrow(props[init_idx]) == 0 && init_idx < length(props)
            init_idx += 1
        end
    end
    trajectories = props[init_idx]
    _start_new_trajectory!(trajectories)

    for candidates in props[2:end]
        # Note: assumes each property table comes from a single observation time!
        nrow(candidates) > 0 && begin
            trajectory_heads = _get_trajectory_heads(
                trajectories, candidates[1, :passtime], maximum_time_step
            )
            candidates_transformed = transform(candidates, candidate_transform_functions...)
            candidates = for floe in eachrow(trajectory_heads)
                candidates_subset = deepcopy(candidates)
                filter_function(floe, candidates_subset)
                nrow(candidates_subset) > 0 && begin
                    candidates_subset[!, :head_uuid] .= floe.uuid
                    candidates_subset[!, :trajectory_uuid] .= floe.trajectory_uuid
                    append!(candidate_pairs, eachrow(candidates_subset))
                end
            end

            # matching function will find best pairs (head_uuid, uuid)
            # and ensure that all pairs are unique
            matched_pairs = DataFrame(candidate_pairs) |> matching_function

            # Get unmatched floes in day 2 (iterations > 2)
            # This should handle the case where there are no trajectory heads available
            matched_uuids = matched_pairs.uuid
            unmatched = filter((f) -> !(f.uuid in matched_uuids), candidates)
            _start_new_trajectory!(unmatched)

            # Attach new matches and unmatched floes to trajectories
            trajectories = vcat(trajectories, matched_pairs, unmatched; cols=:union)
        end
    end
    trajectories = _drop_short_trajectories(trajectories, :trajectory_uuid)
    DataFrames.sort!(trajectories, [:trajectory_uuid, :passtime])
    _add_integer_id!(trajectories, :trajectory_uuid, :ID)
    # Move ID columns to the front
    select!(trajectories, :ID, :trajectory_uuid, :head_uuid, :uuid, :)
    return trajectories
end