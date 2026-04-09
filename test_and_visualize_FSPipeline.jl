using IceFloeTracker
using Images

# setup - running a single case, seeing whether the results are normal
begin
    dataloc = "/Users/dmw/Documents/research/calval_tgrs/data/MODIS_JOG_example_case/"
    falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor")))
    truecolor_files = replace.(falsecolor_files, ("falsecolor" => "truecolor"))
    landmask_file = joinpath(dataloc, "landmask.tiff")
    landmask_img = Gray.(load(landmask_file))
    # fc_file = falsecolor_files[4]
    # tc_file = truecolor_files[4]
    for (tc_file, fc_file) in zip(truecolor_files, falsecolor_files)
        tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
        fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

        segment = FSPipeline.Segment(tile_size_pixels=1200)
        @time begin
            segmentation_results = segment(tc_img, fc_img, landmask_img .> 0)
        end

        save(joinpath(dataloc, "FSPipeline", "binary",
            replace(tc_file, ("truecolor.250m.tiff" => "binary_floes.png")), 
            ), Gray.(labels_map(segmentation_results) .> 0))

        # save colorized image
        cview = view_seg_random(segmentation_results);
        idx = labels_map(segmentation_results) .> 0;
        overlay = deepcopy(tc_img);
        overlay[idx] .= cview[idx];
        save(joinpath(dataloc, "FSPipeline", "colorized",
            replace(tc_file, ("truecolor.250m.tiff" => "colorized_floes.png")), 
            ), overlay)
    end
end

# exports
# save binary image
save(joinpath(dataloc, "FSPipeline", "binary",
    replace(tc_file, ("truecolor.250m.tiff" => "binary_floes_5.png")), 
    ), Gray.(labels_map(test_full) .> 0))
    

function clean_binary_floes(bw_img::AbstractArray{Bool}; min_opening_area=50, strel=FSPipeline.se_disk(3), min_object_size=50) # 
    # erode to separate objects
    out = erode(bw_img, strel)
    hbreak!(out)
    out .= imfill(out, (0, min_object_size)) 
    out .= closing(bw_img)
    # out .= dilate(out, strel)
    out .= .!imfill(.!out, (0, min_opening_area))
    return out
end

tiles = get_tiles(tc_img, 800);
img = Gray.(kmeans_binarization(Gray.(tc_img), fc_img, tiles; cluster_selection_algorithm=IceDetectionBrightnessPeaksMODIS721(0.2, 0.3, 256, 0.01, 3, "union")));
img[1000:2000, 2000:3000]

test_clean = clean_binary_floes(img .> 0);
regions = [(1000:2000, 2000:3000), (2000:3000, 1500:2500)];
mosaicview(mosaicview([img[r...] for r in regions], nrow=1),
mosaicview([tc_img[r...] for r in regions], nrow=1), nrow=2)

# save colorized image
cview = view_seg_random(test_full);
idx = labels_map(test_full) .> 0;
overlay = deepcopy(tc_img);
overlay[idx] .= cview[idx];



binary_overlay = RGB.(Gray.(img));
binary_overlay[idx] .= cview[idx];

save(joinpath(dataloc, "FSPipeline", "colorized",
    replace(tc_file, ("truecolor.250m.tiff" => "colorized_floes_disk5.png")), 
    ), overlay)
save(joinpath(dataloc, "FSPipeline", "colorized",
    replace(tc_file, ("truecolor.250m.tiff" => "colorized_floes_binary_disk5.png")), 
    ), binary_overlay)


Gray.(idx)


mosaicview([view_seg_random(test_full)[r...] for r in regions], nrow=1)

#### Intermediate steps

# Preprocessing

# K-means binarization

# 






p = (diffusion_algorithm=PeronaMalikDiffusion(0.1, 0.1, 5, "exponential"),
    adapthisteq_params=(nbins=256, rblocks=8, cblocks=8, clip=0.01),
    unsharp_mask_params=(smoothing_param=10, intensity=0.5),
    kmeans_params = (k=4, maxiter=50, random_seed=45),
    cluster_selection_algorithm=LopezAcosta2019.IceDetectionLopezAcosta2019())

truecolor_image = float64.(tc_img)
falsecolor_image = float64.(fc_img)
landmask = landmask_img .> 0
coastal_buffer_mask = dilate(landmask_img .> 0, strel_box((51, 51)))

@info "Building cloudmask"
# TODO: Make sure tests aren't over-sensitive to roundoff errors for Float32 vs Float64
cloudmask = create_cloudmask(falsecolor_image)

# 2. Intermediate images
fc_masked = apply_landmask(falsecolor_image, coastal_buffer_mask)

@info "Preprocessing truecolor image"
# nonlinear diffusion
apply_landmask!(truecolor_image, landmask);
sharpened_truecolor_image = nonlinear_diffusion(
    truecolor_image, p.diffusion_algorithm
);

# changed to AdaptiveEqualization directly, could be an issue with channelwise adapthisteq
sharpened_truecolor_image .= IceFloeTracker.Filtering.channelwise_adapthisteq(sharpened_truecolor_image;
    nbins=p.adapthisteq_params.nbins,
    rblocks=p.adapthisteq_params.rblocks,
    cblocks=p.adapthisteq_params.cblocks,
    clip=p.adapthisteq_params.clip
);
sharpened_grayscale_image = unsharp_mask(
    Gray.(sharpened_truecolor_image),
    p.unsharp_mask_params.smoothing_param,
    p.unsharp_mask_params.intensity,
);
apply_landmask!(sharpened_grayscale_image, coastal_buffer_mask)

ice_water_discrim = LopezAcosta2019.discriminate_ice_water(
    sharpened_grayscale_image, fc_masked, coastal_buffer_mask, cloudmask
);

kmeans_result = kmeans_binarization(
            ice_water_discrim,
            fc_masked;
            k=p.kmeans_params.k,
            maxiter=p.kmeans_params.maxiter,
            random_seed=p.kmeans_params.random_seed,
            cluster_selection_algorithm=p.cluster_selection_algorithm
            ) |> LopezAcosta2019.clean_binary_floes

            # check: are there any regions that are nonzero under the cloudmask, since it was applied in discriminate ice water?
apply_cloudmask!(kmeans_result, cloudmask);

# The clean binary floes method has an aggressive fill_holes algorithm. Potentially merging with the
# ice brightness threshold can prevent some of the interstitial water areas from being filled.

# segmentation_B
@info "Segmenting floes part 2/3"
segB = LopezAcosta2019.segmentation_B(sharpened_grayscale_image, cloudmask, kmeans_result);

# Process watershed in parallel using Folds
@info "Building watersheds" # This takes about 5 minutes to run.
@time begin
    watersheds_segB = [
        LopezAcosta2019.watershed_ice_floes(segB.not_ice_bit), LopezAcosta2019.watershed_ice_floes(segB.ice_intersect)
    ]
    watersheds_segB_product = LopezAcosta2019.watershed_product(watersheds_segB...)
end
nothing

# This section takes only 14 seconds to run.
@time begin
segF = LopezAcosta2019.segmentation_F(
    segB.not_ice,
    segB.ice_intersect,
    watersheds_segB_product,
    fc_masked,
    cloudmask,
    coastal_buffer_mask,
)
end

# When I run it in here, the results look good.
# When I run it in the other script, it looks bad. Not sure how that can be!

