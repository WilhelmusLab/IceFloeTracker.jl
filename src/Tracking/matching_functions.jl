

"""
    MinimumWeightMatchingFunction(columns=[:scaled_distance, :relative_error_area, ...])
    MinimumWeightedMatchingFunction(candidate_pairs)

Function to identify a best matching between pairs of ice floes in the DataFrame `candidate_pairs`. The
`columns` variable is instantiated by the first functor call and is used to select a list of columns in `candidate_pairs`
to sum. The result of the sum is the weight assigned to each pairing. Then, a best set of unique pairs is found by
carrying out two grouped minimizations: first grouping by the first floe, identified by the `head_uuid` column, and
finding the floe with the smallest weight, then grouping by the second floe, identified by the `uuid` column, and again 
finding the floe with the smallest weight. Finally, only pairs that exist in both the forward and backward grouped minizations
are identified as likely true matches.
"""
@kwdef struct MinimumWeightMatchingFunction <: AbstractFloeFilterFunction
    columns=[:scaled_distance, :relative_error_area, :relative_error_convex_area, 
                    :relative_error_major_axis_length, :relative_error_minor_axis_length,
                    :psi_s_correlation_score, :scaled_shape_difference]
end

function (f::MinimumWeightMatchingFunction)(candidate_pairs::DataFrame);
     
    # Potential future updates: replace sum with a weighted
    candidate_pairs[!, :w] = sum.(eachrow(candidate_pairs[:, f.columns]))

    # Forward: f -> {g}, find minimum dx over set {g}
    matches_fwd = combine(sdf -> sdf[argmin(sdf.w), :], groupby(candidate_pairs, :head_uuid));
    matches_fwd = combine(sdf -> sdf[argmin(sdf.w), :], groupby(candidate_pairs, :uuid));

    # Backward: {f} <- g, find minimum dx over {f}
    matches_bwd = combine(sdf -> sdf[argmin(sdf.w), :], groupby(candidate_pairs, :uuid));
    matches_bwd = combine(sdf -> sdf[argmin(sdf.w), :], groupby(candidate_pairs, :head_uuid));
    
    return innerjoin(matches_fwd[:, [:head_uuid, :uuid]], matches_bwd, on = [:head_uuid, :uuid]);
end


# """
#     find_floe_matches(
#     tracked,
#     candidate_props,
#     condition_thresholds,
#     mc_thresholds
# )

# Find matches for floes in `tracked` from floes in  `candidate_props`.

# # Arguments
# - `tracked`: dataframe containing floe trajectories.
# - `candidate_props`: dataframe containing floe candidate properties.
# - `candidate_filter_settings`: thresholds for deciding whether to match floe `i` from tracked to floe j from `candidate_props`
# - `candidate_matching_settings`: thresholds for area mismatch and psi-s shape correlation
# """
# function find_floe_matches(
#     tracked::T, candidate_props::T, candidate_filter_settings, candidate_matching_settings
# ) where {T<:AbstractDataFrame}
#     matches = []
#     for floe1 in eachrow(tracked), floe2 in eachrow(candidate_props)
#         Δt = floe2.passtime - floe1.passtime
#         ratios, conditions, dist = compute_ratios_conditions(
#             floe1, floe2, Δt, candidate_filter_settings
#         )

#         if callmatchcorr(conditions)
#             (area_mismatch, corr) = matchcorr(
#                 floe1.mask, floe2.mask, Δt; candidate_matching_settings.comp...
#             )
#             if isfloegoodmatch(conditions, candidate_matching_settings.goodness, area_mismatch, corr)
#                 @debug "** Found a good match for ", floe1.uuid, "<=", floe2.uuid
#                 push!(
#                     matches,
#                     (;
#                         Δt,
#                         measures=(; ratios..., area_mismatch, corr),
#                         conditions,
#                         dist,
#                         floe1,
#                         floe2,
#                     ),
#                 )
#             end
#         end
#     end

#     remaining_matches_df = DataFrame(matches)
#     best_matches = []

#     for floe2 in eachrow(
#         sort(candidate_props, :area; rev=true),  # prioritize larger floes
#     )
#         matches_involving_floe2_df = filter((r) -> r.floe2 == floe2, remaining_matches_df)
#         nrow(matches_involving_floe2_df) == 0 && continue
#         best_match = matches_involving_floe2_df[1, :]
#         measures_df = DataFrame(matches_involving_floe2_df.measures)
#         best_match_idx = getidxmostminimumeverything(measures_df) # 
#         best_match = matches_involving_floe2_df[best_match_idx, :]
#         push!(
#             best_matches,
#             (;
#                 best_match.floe2...,
#                 area_mismatch=best_match.measures.area_mismatch,
#                 corr=best_match.measures.corr,
#                 trajectory_uuid=best_match.floe1.trajectory_uuid,
#                 head_uuid=best_match.floe1.uuid,
#             ),
#         )
#         # Filter out from the remaining matches any cases involving the matched floe1
#         remaining_matches_df = filter(
#             (r) -> !(r.floe1 === best_match.floe1), remaining_matches_df
#         )
#         # Filter out from the remaining matches any cases involving the matched floe2
#         remaining_matches_df = filter(
#             (r) -> !(r.floe2 === best_match.floe2), remaining_matches_df
#         )
#     end

#     best_matches_df = similar(tracked, 0)
#     append!(best_matches_df, best_matches; promote=true)
#     return best_matches_df
# end
