begin
    HOME = "." # path to the root of the project two levels up

    # Activate the environment
    # using Pkg
    # Pkg.activate(HOME)
    # Pkg.precompile()

    using IceFloeTracker: pairfloes, deserialize, PaddedView, float64, mosaicview, Gray, sort_floes_by_area!, add_passtimes!, addfloemasks!, addψs!, getpropsday1day2, MatchedPairs, makeemptydffrom, compute_ratios_conditions, callmatchcorr, isfloegoodmatch, appendrows!, getidxmostminimumeverything, isnotnan, getbestmatchdata, addmatch!, resolvecollisions!, deletematched!, update!, Tracked, get_trajectory_heads, get_unmatched, _pairfloes, _swap_last_values!, get_dt, find_floe_matches, consolidate_matched_pairs
    using DataFrames
    using Dates
    using Random
    using StatsBase
    using Test
    imshow(x) = Gray.(x)
end


begin

    begin # Set thresholds
        t1 = (dt=(30.0, 100.0, 1300.0), dist=(200, 250, 300))
        t2 = (
            area=1200,
            arearatio=0.28,
            majaxisratio=0.10,
            minaxisratio=0.12,
            convexarearatio=0.14,
        )
        t3 = (
            area=10_000,
            arearatio=0.18,
            majaxisratio=0.1,
            minaxisratio=0.15,
            convexarearatio=0.2,
        )
        condition_thresholds = (t1, t2, t3)
        mc_thresholds = (
            goodness=(area3=0.18, area2=0.236, corr=0.68), comp=(mxrot=10, sz=16)
        )
    end


    begin # Load data
        pth = joinpath(HOME, "test", "test_inputs", "tracker")
        _floedata = deserialize(joinpath(pth, "tracker_test_data.dat"))
        _passtimes = deserialize(joinpath(pth, "passtimes.dat"))
        latlonimgpth = joinpath(HOME, "test", "test_inputs", "NE_Greenland_truecolor.2020162.aqua.250m.tiff")
        _props, _imgs = deepcopy.([_floedata.props, _floedata.imgs])

        # This order is important: masks, uuids, passtimes, ψs
        IceFloeTracker.addfloemasks!(_props, _imgs)
        IceFloeTracker.adduuid!(_props)
        IceFloeTracker.add_passtimes!(_props, _passtimes)
        IceFloeTracker.addψs!(_props)
    end

    begin # Filter out floes with area less than `floe_area_threshold` pixels
        floe_area_threshold = 400
        for (i, prop) in enumerate(_props)
            _props[i] = prop[prop[:, :area].>=floe_area_threshold, :] # 500 working good
            sort!(_props[i], :area, rev=true)
        end
    end

    # Prep data for _pairfloes. TODO: use as test case
    begin
        # Unmatched floe in day 1, unmatched floe in day 2, and matches for every floe starting in day 3
        props_test_case2 = deepcopy(_props)
        delete!(props_test_case2[1], 1)
        delete!(props_test_case2[2], 5)
        segmented_imgs = _imgs
        passtimes = _passtimes
        # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2
    end

    begin # 0th iteration: pair floes in day 1 and day 2 and add unmatched floes to _pairs
        matched_pairs = find_floe_matches(props_test_case2[1], props_test_case2[2], condition_thresholds, mc_thresholds)

        # Get unmatched floes from day 1
        unmatched1 = get_unmatched(props_test_case2[1], matched_pairs.props1)
        unmatched2 = get_unmatched(props_test_case2[2], matched_pairs.props2)
        unmatched = vcat(unmatched1, unmatched2)
        consolidated_matched_pairs = consolidate_matched_pairs(matched_pairs)
        consolidated_matched_pairs[:, [:uuid, :passtime, :area_mismatch, :corr]]
    end


    begin # Update _pairs with unmatched floes
        _pairs = vcat(consolidated_matched_pairs, unmatched[:, names(consolidated_matched_pairs)])
        _pairs[:, [:uuid, :passtime, :area_mismatch, :corr]]
    end

    begin # Start 2:end iterations
        trajectory_heads = get_trajectory_heads(_pairs)
        new_pairs = IceFloeTracker.find_floe_matches(trajectory_heads, props_test_case2[3], condition_thresholds, mc_thresholds)
        # Get unmatched floes in day 2 (iterations > 2)
        unmatched2 = get_unmatched(props_test_case2[3], new_pairs.props2)
        @assert isempty(unmatched2)

        new_pairs = IceFloeTracker.get_matches(new_pairs)
    end

    _pairs = vcat(_pairs, new_pairs, unmatched2)
    DataFrames.sort!(_pairs, [:uuid, :passtime])
    _swap_last_values!(_pairs)
    IceFloeTracker.reset_id!(_pairs)
    _pairs[:, [:uuid, :passtime, :area_mismatch, :corr]]
end

# Test adding an unmatched floe in props2
begin
    fake_props2 = deepcopy(props[3])
    newrow = deepcopy(fake_props2[end, :])
    newrow.area = 5
    uuidfake = "zfakefloe123"
    newrow.uuid = uuidfake
    push!(fake_props2, newrow)
    unmatched_df = get_unmatched(fake_props2, foo.props2)
    @assert unmatched_df[1, :uuid] == uuidfake
    @assert "corr" in names(unmatched_df)
end

# Repeat

# At the end make :uuid => :ID where ID is a unique identifier for each floe from 1 to length(unique(:uuid))