# Run time and error rate with the morph floes vs using watershed + final morph cleanupusing DataFrames
# Before re-running, investigate the errors in the first tests.
# -- Why aren't floes being split?
# -- Why are very large floes retained?
# -- Would it work to improve the distance split function?
# -- Would it work to improve the morph-split function?
# -- How are the linear floe separation artifacts showing up?

using Images
using IceFloeTracker

dataset = Watkins2026Dataset(; ref="v0.1")
cases = filter(f -> f.visible_floes == "yes" && f.region == "baffin_bay", dataset)
cases = filter(f -> f.case_number == 128, dataset)
segment = Watkins2026.Segment(;tile_size_pixels=200)

segment = LopezAcosta2019.Segment()
# segment = LopezAcosta2019Tiling.Segment()
# Morphed split failed for Case 63 - caresian index error!
# Case 61 - highlights error with partial cloud cover. Could flag floes that are bordered by clouds, and have
# special case for tracking for partial obscuration.
# Case 128 doesn't work well - near ice edge, less distinct ice edges. I think the boundary gradient would be enough to remove them, but a better
# preprocessing routine could make the difference.
@time begin
    segmented_images = segment.(modis_truecolor.(cases), modis_falsecolor.(cases), modis_landmask.(cases));
end

# Overlay outlines and show with original image
# Takes 30 seconds for the first pair of images, then running other ones is faster.
viz = []
for (case, segments) in zip(cases, segmented_images)
    outlines = validated_binary_floes(case)
    outlines .= outlines .- erode(outlines)

    cview_overlay = RGB.(modis_truecolor(case))
    cview_final = view_seg_random(segments)
    cview_overlay[labels_map(segments) .> 0] .= cview_final[labels_map(segments) .> 0]
    cview_final[labels_map(segments) .== 0] .= RGB(0,0,0)
    cview_final[outlines .> 0] .= RGB(1,1,1)
    push!(viz, mosaicview(modis_truecolor(case), cview_overlay, cview_final, nrow=1))
end
mosaicview([v for v in viz], nrow=2)

#### Example Cases from April 2014 ####
# Something broke on the LopezAcosta2019 case, so I need to dig into that. Previously it had worked well.
# Re-running the Watkins2026 method, using the updated floe-split routine. Started at 11:50 ish.
begin
    # segment = Watkins2026.Segment(;tile_size_pixels=1000)
    dataloc = "/Users/dmw/Documents/research/calval_tgrs/data/MODIS_JOG_example_case/"
    falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor")))
    truecolor_files = replace.(falsecolor_files, ("falsecolor" => "truecolor"))
    landmask_file = joinpath(dataloc, "landmask.tiff")
    landmask_img = Gray.(load(landmask_file));
end

for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    @info "Processing image " * tc_file
    @time begin
        segment_result = segment(tc_img, fc_img, landmask_img)
    end

    # save binary image
    save(joinpath(dataloc, "LopezAcosta2019", "binary",
        replace(tc_file, ("truecolor.250m.tiff" => "binary_floes_.png")), 
        ), Gray.(labels_map(segment_result) .> 0))
        
    # save colorized image
    cview = view_seg_random(segment_result)
    idx = binarize_segments(segment_result) .> 0
    tc_img[idx] .= cview[idx]
    
    save(joinpath(dataloc, "LopezAcosta2019", "colorized",
        replace(tc_file, ("truecolor.250m.tiff" => "colorized_floes.png")), 
        ), tc_img)
end


###### Potential speedups ######
# Lazy evaluation of solidity (only compute after checking area and circularity)
# Check size of thin connections using local minima of the distance transform
bw = labels_map(segmented_images[1]) .== 0;
Gray.(bw)
d = .- distance_transform(feature_transform(bw));
0.5 .* Gray.(bw) .+ Gray.(dilate(local_minima(d), strel_diamond((5,5))))

indices = component_indices(labels_map(segmented_images[1]))
boxes = component_boxes(labels_map(segmented_images[1]))
effective_radii = Dict(r => sqrt(area[r]/pi) for r in keys(area));
candidates = Dict()
for r in keys(indices)
    if r > 0
        dist_list = d[indices[r]]
        minima_list = (local_minima(d) .> 0)[indices[r]]
        push!(candidates, r => abs.(unique(round.(Int64, dist_list[minima_list]))))
    end
end
masks = component_floes(labels_map(segmented_images[1]));
component_convex_areas(Int64.(masks[1]))
indexmap = labels_map(segmented_images[1])
areas = component_lengths(indexmap)
perimeters = component_perimeters(indexmap)
masks = component_floes(indexmap)
# convex_areas = component_convex_areas(indexmap)
labels = filter(r -> r != 0, intersect(keys(areas), keys(perimeters)))

circularities = Dict(r => 4 * pi * areas[r] / perimeters[r]^2 for r in labels if r != 0)


##### Checking iterations of the morph split function #####
# - Update doesn't work as well as it did

max_depth=5
min_area=100
min_circularity=0.6
min_solidity=0.85
max_iter=10


# apply opening with radius r = 1:max_depth until either the
# shape is broken into parts or you reach max depth. If max depth
# is reached and nothing changes, return the original mask.
# To do: figure out how to keep going until the large objects have been split enough times
function se_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    
    # _generate_se!(se)
    return sqrt.(se) .<= r
end

function morph_split(mask)
    _max_depth = round(Int, sqrt(sum(mask)/pi))
    _mask = deepcopy(mask)
    for r in 1:_max_depth
        new_labels = opening(_mask, se_disk(r)) |> label_components
        length(unique(new_labels)) > 2 && return new_labels .> 0
    end
    return nothing
end

# Get the list of candidates for splitting by checking area, circularity, and solidity.
# Evaluates solidity lazily since that operation is slow.
function get_candidate_labels(indexmap, min_area, min_circularity, min_solidity)
    areas = component_lengths(indexmap)
    perimeters = component_perimeters(indexmap)
    masks = component_floes(indexmap)
    # convex_areas = component_convex_areas(indexmap)
    labels = filter(r -> r != 0, intersect(keys(areas), keys(perimeters)))
    
    circularities = Dict(r => 4 * pi * areas[r] / perimeters[r]^2 for r in labels if r != 0)
    
    # Only evaluate convex area if the area and circularity pass the thresholds for splitting
    split_labels = Int64[]
    for r in labels
        (r != 0) && (areas[r] > 1.5 * min_area) && (circularities[r] < min_circularity) && begin
            ca = component_convex_areas(Int64.(masks[r]))
            s = areas[r] / ca[1]
            s < min_solidity && push!(split_labels, r)
        end
    end
    return split_labels
end

begin # Initialize
    indexmap = label_components(labels_map(segmented_images[1]))
    out = deepcopy(indexmap)

    max_depth = 10
    label_offset = maximum(indexmap)
    boxes = component_boxes(indexmap)
    indices = component_indices(indexmap)
    masks = component_floes(indexmap)

    n_regions = length(unique(out))
    n_updated = 0
    split_labels = get_candidate_labels(out, min_area, min_circularity, min_solidity)
    count = 0
    out_img = 0.5 .* Gray.(out .> 0)
    # out_img[isboundary(out) .> 0] .= Gray(0)
    for r in split_labels
        out_img[indices[r]] .= Gray(1)
    end
    intermediate_images = [out_img]
end

while (n_regions != n_updated) && (count < max_iter)
    n_regions = length(unique(out))
    for r in split_labels
        update_floe = morph_split(masks[r])
        !isnothing(update_floe) && begin
            update_labels = label_components(update_floe)
            update_areas = component_lengths(update_labels)
            
            for rnew in keys(update_areas)
                # imfill works for fully separated components, but will miss if components are touching.                
                update_areas[rnew] < min_area && (update_labels[update_labels .== rnew] .= 0)
            end
            out[indices[r]] .= 0 # Remove original floe
            out[boxes[r]] .+= update_labels # New floe has at most equal area to the original, so there shouldn't be overlap
        end
    end

    out .= label_components(out)
    count += 1
    n_updated = length(unique(out))
    println(count, " ", n_regions, " ", n_updated)
    n_regions != n_updated && begin 
        
        # Select candidate labels which were successfully split in the morph split step.
        split_labels = [r for r in get_candidate_labels(out, min_area, min_circularity, min_solidity) 
                            if length(unique(out[indices[r]])) > 1]

        boxes = component_boxes(out)
        indices = component_indices(out)
        masks = component_floes(out)

        # overlay new list of labels to split
        out_img = 0.5 .* Gray.(out .> 0)
        for r in split_labels
            out_img[indices[r]] .= Gray(1)
        end
        push!(intermediate_images, deepcopy(out_img))
    end
end


@info "Morph split completed in "*string(count)*" iterations"
final_indexmap = label_components(out)

mosaicview(intermediate_images..., ncol=3, rowmajor=true)

0.5 .* Gray.(intermediate_images[1] .> 0) .+ Gray.(intermediate_images[end] .> 0)


# Finding reason for intermediate image fail
tc_file = truecolor_files[2]
fc_file = falsecolor_files[2]

img = Gray.(load(joinpath(dataloc, "Watkins2026", "binary",
        replace(tc_file, ("truecolor.250m.tiff" => "binary_floes.png")), 
        )));
indexmap = label_components(img);
tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)));
fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)));
land_mask = landmask_img .> 0
cloud_mask = Watkins2025CloudMask()(fc_img)
cloud_mask .= apply_landmask(cloud_mask, land_mask)
Gray.(cloud_mask)

prelim_ice_mask = Watkins2026.ice_water_mask(tc_img, cloud_mask, land_mask);

Gray.(prelim_ice_mask) .+ 0.5 .* Gray.(cloud_mask)


bright_ice = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01
    )(apply_landmask(fc_img, landmask_img .> 0))
nothing