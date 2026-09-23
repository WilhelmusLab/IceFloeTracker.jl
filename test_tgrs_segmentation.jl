using IceFloeTracker
using Images
using DataFrames
using Dates

dataloc = "workflow/results/filetype"
case = "baffin_bay_006.250m.2022-05-30.aqua.tiff"

tc_imgs = load.([joinpath(dataloc, "truecolor", case), joinpath(dataloc, "truecolor", replace(case, "aqua" => "terra"))])
fc_imgs = load.([joinpath(dataloc, "falsecolor", case), joinpath(dataloc, "falsecolor", replace(case, "aqua" => "terra"))])
lm_imgs = load.([joinpath(dataloc, "landmask", case), joinpath(dataloc, "landmask", replace(case, "aqua" => "terra"))]);
land_masks = lm_imgs .|> r -> Gray.(r) .> 0;

mosaicview(tc_imgs[1], fc_imgs[1], nrow=1) 

preprocess = TGRS2026.Preprocess();
preproc_gray = preprocess.(tc_imgs, land_masks);

classify = TGRS2026.Classify();
classified = classify.(fc_imgs, land_masks);

mosaicview(TGRS2026.colorize_classification.(classified), nrow=1)
@time begin
test = TGRS2026.Segment().(tc_imgs, fc_imgs, land_masks);
end

cviews = view_seg_random.(test);
overlay_outlines = TGRS2026.colorize_classification.(classified)
for idx in 1:2
    bdry = isboundary(labels_map(test[idx]))
    overlay_outlines[idx][bdry .> 0] .= cviews[idx][bdry .> 0]
end
mosaicview(overlay_outlines..., nrow=1)


overlay_tc = copy(tc_imgs)
for idx in 1:2
    bdry = isboundary(labels_map(test[idx]))
    overlay_tc[idx][bdry .> 0] .= cviews[idx][bdry .> 0]
end
mosaicview(overlay_tc..., nrow=1)
save("../../test_006_overlay_2.png", mosaicview(overlay_tc..., nrow=1))

@time begin
    fmt = DateFormat("y-m-dTH:M:S")
    pass_times = [DateTime("2022-05-30T15:28:46", fmt),
                  DateTime("2013-03-08T17:56:22", fmt)]
    floe_tracker = TGRS2026.Track(minimum_area=100);
    tracked_floes = floe_tracker(test, pass_times);
end

# TODO: set up filter to flag floes which are not tracked and above minimum tracking area


