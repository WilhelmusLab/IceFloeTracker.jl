# This script is to test out an updated, streamlined LA2019 approach I'm calling FSPipeline
# for now. It's the LA2019 algorithm with tiling. I've replaced the watershed function with a
# simpler one, with the marker made with the distance transform cut off at 4 pixel depth and then
# eroded.

# The algorithm is much faster, but at some point I'm getting the same issue I had earlier
# where the whole area gets turned into a giant floe. It's not the cluster selection algorithm.
# The k-means result looks wrong.


using IceFloeTracker
using Images

# setup
begin
    dataloc = "/Users/dmw/Documents/research/calval_tgrs/data/MODIS_JOG_example_case/"
    falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor")))
    truecolor_files = replace.(falsecolor_files, ("falsecolor" => "truecolor"))
    landmask_file = joinpath(dataloc, "landmask.tiff")
    landmask_img = Gray.(load(landmask_file))
    # for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)

    fc_file = falsecolor_files[1]
    tc_file = truecolor_files[1]
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    # segment = LopezAcosta2019.Segment()
    segment = FSPipeline.Segment(
        cloud_mask_algorithm=LopezAcostaCloudMask(),
        coastal_buffer_structuring_element=strel_box((5,5))
    )
end

# 
# zoom = (1000:2500, 1000:2500)
# Already looks different (and bad) after the k-means step
# The discriminate ice water function seems to rely on the
# image being strongly equalized (0.99 didn't work, 0.01 seemed to
# actually be fine.)

segment = FSPipeline.Segment(
                cloud_mask_algorithm=LopezAcostaCloudMask(),
                coastal_buffer_structuring_element=strel_box((5,5)),
                adapthisteq_params=(nbins=256, rblocks=8, cblocks=8, clip=0.2),
                )

for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    @time begin
        segmented_image = segment(tc_img, fc_img, landmask_img)
    end

    save(joinpath(dataloc, "LopezAcosta2019", "binary",
        replace(tc_file, ("truecolor.250m.tiff" => "binary_floes_tiled_strong_equalization.png")), 
        ), Gray.(segmented_image.image_indexmap .> 0))
end

# testing: exporting just the kmeans result.
Gray.(test)
# Gray.(test.image_indexmap .> 0)
# check what the resolution was on the edges in the watershed boundary

p = (diffusion_algorithm=PeronaMalikDiffusion(0.1, 0.1, 5, "exponential"),
    adapthisteq_params=(nbins=256, rblocks=8, cblocks=8, clip=0.2),
    unsharp_mask_params=(smoothing_param=10, intensity=0.5),
    kmeans_params = (k=4, maxiter=50, random_seed=45),
    tile_size_pixels = 1000,
    cloud_mask_algorithm = Watkins2025CloudMask(),
    min_tile_ice_pixel_count = 1000,
    coastal_buffer_structuring_element = strel_box((51,51)),
    cluster_selection_algorithm=LopezAcosta2019.IceDetectionLopezAcosta2019())


for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    @time begin
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))
    end    
    break
end

begin
    @info "Setting up initial masks"
    begin
        # TODO: Make sure tests aren't over-sensitive to roundoff errors for Float32 vs Float64
        _lm_temp = create_landmask(landmask_img, p.coastal_buffer_structuring_element)
        land_mask = _lm_temp.non_dilated
        coastal_buffer = _lm_temp.dilated
        fc_masked = apply_landmask(fc_img, coastal_buffer)
        cloud_mask = create_cloudmask(fc_masked) # this way the cloud mask is already land masked
        apply_cloudmask!(fc_masked, cloud_mask)
    

        # Get tiles
        n, m = size(tc_img)
        tile_size_pixels = p.tile_size_pixels
        tile_size_pixels > maximum([n, m]) && begin
            @warn "Tile size too large, defaulting to image size"
            tile_size_pixels = minimum([n, m])
        end
    end

    tiles = get_tiles(tc_img, tile_size_pixels);
    cloud_mask = p.cloud_mask_algorithm(fc_img);

    # very basic sanity check to find sea ice pixels
    prelim_ice_mask = (red.(tc_img) .> 0.3) .&& .!cloud_mask .&& .!coastal_buffer
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > p.min_tile_ice_pixel_count, # check if the min is larger than the tile size and adjust
        tiles);

    @info "Preprocessing truecolor image"
    # nonlinear diffusion
    apply_landmask!(tc_img, land_mask);
    preproc_img = float64.(deepcopy(tc_img));

    for tile in filtered_tiles
        preproc_img[tile...] .= nonlinear_diffusion(
            preproc_img[tile...], p.diffusion_algorithm
        )
    end

    # changed to AdaptiveEqualization directly, could be an issue with channelwise adapthisteq
    # Full image, but we'll zero out the masks again
    preproc_img .= IceFloeTracker.Filtering.channelwise_adapthisteq(preproc_img;
        nbins=p.adapthisteq_params.nbins,
        rblocks=p.adapthisteq_params.rblocks,
        cblocks=p.adapthisteq_params.cblocks,
        clip=p.adapthisteq_params.clip
    );

    preproc_gray = Gray.(preproc_img)
    for tile in filtered_tiles
        preproc_gray[tile...] .= unsharp_mask(
        preproc_gray[tile...],
        p.unsharp_mask_params.smoothing_param,
        p.unsharp_mask_params.intensity)
    end

    apply_landmask!(preproc_gray, coastal_buffer)
    apply_cloudmask!(preproc_gray, cloud_mask);
    nothing
end

# Preprocessing seems fine - maybe a little too intense with the sharpening, but not crazy.
Gray.(preproc_gray)
save("/Users/dmw/Downloads/preproc_strong_eq.png", preproc_gray)

test_adapt_bin = tiled_adaptive_binarization(preproc_gray, filtered_tiles; 
minimum_window_size=400, threshold_percentage=0, minimum_brightness=0.4);

save("/Users/dmw/Downloads/preproc_strong_eq_binarized.png", Gray.(test_adapt_bin))


begin
    ice_water_discrim = zeros(size(preproc_gray))
    for tile in filtered_tiles
        ice_water_discrim[tile...] .= LopezAcosta2019.discriminate_ice_water(
            preproc_gray[tile...], fc_masked[tile...], coastal_buffer[tile...], cloud_mask[tile...]
    );
    end
end
alt_iwd = deepcopy(test)
mosaicview(Gray.(ice_water_discrim[zoom...]), Gray.(alt_iwd[zoom...]), nrow=1)

begin
    kmeans_result = kmeans_binarization(
                Gray.(alt_iwd), # errors if not Gray type
                fc_masked, 
                filtered_tiles;
                k=p.kmeans_params.k,
                maxiter=p.kmeans_params.maxiter,
                random_seed=p.kmeans_params.random_seed,
                cluster_selection_algorithm=p.cluster_selection_algorithm
                ) |> LopezAcosta2019.clean_binary_floes   
    apply_landmask!(kmeans_result, coastal_buffer);
    nothing
end
zoom = (1000:2000, 1000:3000)
# fc_masked[zoom...]
mosaicview(Gray.(kmeans_result)[zoom...],
    Gray.(ice_water_discrim)[zoom...])

begin
    # The clean binary floes method has an aggressive fill_holes algorithm. Potentially merging with the
    # ice brightness threshold can prevent some of the interstitial water areas from being filled.

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = LopezAcosta2019.segmentation_B(preproc_gray, cloud_mask, kmeans_result);

    # # Process watershed in parallel using Folds
    # @info "Building watersheds" # This takes about 5 minutes to run in the original version.
    # @time begin
    #     watersheds_segB = [
    #         LopezAcosta2019.watershed_ice_floes(segB.not_ice_bit), LopezAcosta2019.watershed_ice_floes(segB.ice_intersect)
    #     ]
    #     watersheds_segB_product = LopezAcosta2019.watershed_product(watersheds_segB...)
    # end
    # nothing

    # instead, use the watershed transform from the new algorithm
    @info "tiled watershed transform on the intersect"
    @time begin
    w_merged = Watkins2026.watershed_transform(
        segB.ice_intersect,
        segB.not_ice,
        filtered_tiles;
        strel=strel_diamond((5,5)), # input param
        dist_threshold=4 # input param
    )
    end
    @info "tiled watershed on the not ice bit, whatever that one is"
    @time begin
    w_other = Watkins2026.watershed_transform(
        segB.not_ice_bit,
        segB.not_ice,
        filtered_tiles;
        strel=strel_diamond((5,5)),
        dist_threshold=4
    )
    end

    watersheds_segB_product = falses(size(tc_img))
    for tile in filtered_tiles
        watersheds_segB_product[tile...] .= (isboundary(labels_map(w_merged)[tile...]) .* isboundary(labels_map(w_other)[tile...])) .> 0
    end

    # This section takes only 14 seconds to run.
    @time begin
    segF = LopezAcosta2019.segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB_product,
        fc_masked,
        cloud_mask,
        coastal_buffer
    )
    end

    save(joinpath(dataloc, "LopezAcosta2019", "binary",
    replace(tc_file, ("truecolor.250m.tiff" => "binary_floes_tiled.png")), 
    ), Gray.(segF))

end


