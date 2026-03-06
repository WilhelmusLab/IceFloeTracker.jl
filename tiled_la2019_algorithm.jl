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
    # dataloc = "/Users/dmw/Documents/research/calval_tgrs/data/MODIS_JOG_example_case/"
    # dataloc = "/Volumes/Research/ENG_Wilhelmus_Shared/group/IFT_greenland_sea_floe_traje ctories/greenland_sea_ift_dataset/"
    dataloc = "/Users/dwatkin2/Documents/research/manuscripts/greenland_floe_scale_dataset/ift_fram_strait_test_cases/data/modis/"
    # Year is needed for the greenland sea dataset on the share
    # year = string(2003)
    # falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor", year)))
    falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor")))

    case_number = "009" # Next: 011
    falsecolor_files = filter(f -> occursin(case_number, f), falsecolor_files)

    truecolor_files = replace.(falsecolor_files, ("falsecolor" => "truecolor"))
    landmask_file = joinpath(dataloc, "landmask.tiff")
    landmask_img = Gray.(load(landmask_file)) 
    nothing
end


# This gives us 120 images for testing.
# zoom = (1000:2500, 1000:2500)
# Already looks different (and bad) after the k-means step
# The discriminate ice water function seems to rely on the
# image being strongly equalized (0.99 didn't work, 0.01 seemed to
# actually be fine.)

segment = FSPipeline.Segment(
                cloud_mask_algorithm=LopezAcostaCloudMask(),
                coastal_buffer_structuring_element=strel_box((5,5)),
                tile_size_pixels=500,
                adapthisteq_params=(nbins=256, rblocks=11, cblocks=6, clip=20),
                )

# TODO: Add callback function to return the binarized image prior to floe-splitting.
# TODO: Add callback function to return the preprocessed sharpened image and the ice-water-discrimination image

# saveloc = "/Users/dwatkin2/Documents/research/manuscripts/greenland_floe_scale_dataset/IFT_greenland_sea_dataset/data/FSPipeline_results/"
saveloc = "/Users/dwatkin2/Documents/research/manuscripts/greenland_floe_scale_dataset/IFT_greenland_sea_dataset/data/FSPipeline_results/binary_lacm_clip20_test_cases/"

for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file))) # add year if pulling from share
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    @time begin
        segmented_image = segment(tc_img, fc_img, landmask_img)
    end
    save(joinpath(saveloc, 
        replace(tc_file, ("truecolor.tiff" => "binary_floes.png")), 
        ), Gray.(segmented_image.image_indexmap .> 0))
end

segment = LopezAcosta2019.Segment()
for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file))) # add year if pulling from share
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    @time begin
        segmented_image = segment(tc_img, fc_img, landmask_img)
    end
    save(joinpath(saveloc, "LA2019_julia",
        replace(tc_file, ("truecolor.tiff" => "binary_floes.png")), 
        ), Gray.(segmented_image.image_indexmap .> 0))
end


#### In this section: running section by section so I can view the output ####
# testing: exporting just the kmeans result.
# Gray.(test.image_indexmap .> 0)
# check what the resolution was on the edges in the watershed boundary

p = (diffusion_algorithm=PeronaMalikDiffusion(0.1, 0.1, 5, "exponential"),
    adapthisteq_params=(nbins=256, rblocks=8, cblocks=8, clip=0.2),
    unsharp_mask_params=(smoothing_param=10, intensity=0.5),
    kmeans_params = (k=4, maxiter=50, random_seed=45),
    tile_size_pixels = 500,
    cloud_mask_algorithm = Watkins2025CloudMask(),
    min_tile_ice_pixel_count = 1000, # or use rblocks / cblocks
    coastal_buffer_structuring_element = strel_box((51,51)),
    cluster_selection_algorithm=IceDetectionBrightnessPeaksMODIS721(0.1, 0.3, 64, 0.01, 3, "union")
    )


for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    @time begin
    global tc_img = RGB.(load(joinpath(dataloc, "truecolor", year, tc_file)))
    global fc_img = RGB.(load(joinpath(dataloc, "falsecolor", year, fc_file)))
    end    
    break
end

mosaicview(tc_img, fc_img, nrow=1)

@time begin
    init_segment_results = segment(tc_img, fc_img, landmask_img)
end
begin
    # coastal intersection, update buffer
    coastal_buffer = dilate(landmask_img .> 0, strel_box((25,25)))
    binary_segments = labels_map(init_segment_results) .> 0
    potential_landfast = mreconstruct(dilate, coastal_buffer, binary_segments, strel_box((3,3)))
end

mosaicview(
    [Gray.(m) for m in [coastal_buffer, binary_segments, potential_landfast]], nrow=1)

begin
updated_mask = (landmask_img .> 0) .|| potential_landfast
post_segment_results = segment(tc_img, fc_img, updated_mask)
updated_binary_segments = labels_map(post_segment_results) .> 0
end
mosaicview(
    [Gray.(m) for m in [binary_segments, updated_binary_segments]], nrow=1)

save("/Users/dwatkin2/Downloads/original_binary_floes_500px.png", Gray.(binary_segments))
save("/Users/dwatkin2/Downloads/updated_mask_binary_floes_500px.png", Gray.(updated_binary_segments))
save("/Users/dwatkin2/Downloads/tc_img.png", tc_img)
save("/Users/dwatkin2/Downloads/band7.png", Gray.(red.(fc_img)))
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
    # preproc_img .= IceFloeTracker.Filtering.channelwise_adapthisteq(preproc_img;
    #     nbins=p.adapthisteq_params.nbins,
    #     rblocks=p.adapthisteq_params.rblocks,
    #     cblocks=p.adapthisteq_params.cblocks,
    #     clip=p.adapthisteq_params.clip
    # );
    @time begin
        adjust_histogram!(preproc_img,
            ContrastLimitedAdaptiveHistogramEqualization(
                nbins=p.adapthisteq_params.nbins,
                rblocks=p.adapthisteq_params.rblocks,
                cblocks=p.adapthisteq_params.cblocks,
                clip=p.adapthisteq_params.clip)
        )
    end


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


begin
    ice_water_discrim = zeros(size(preproc_gray))
    for tile in filtered_tiles
        ice_water_discrim[tile...] .= LopezAcosta2019.discriminate_ice_water(
            preproc_gray[tile...], fc_masked[tile...], coastal_buffer[tile...], cloud_mask[tile...]
    );
    end
end

# There's an issue
ice_water_discrim[coastal_buffer .> 0] .= Gray(1)
ice_water_discrim[cloud_mask .> 0] .= Gray(1)

Gray.(ice_water_discrim[zoom...])

begin
    kmeans_result = kmeans_binarization(
                Gray.(ice_water_discrim), # errors if not Gray type
                fc_masked, 
                filtered_tiles;
                k=p.kmeans_params.k,
                maxiter=p.kmeans_params.maxiter,
                random_seed=p.kmeans_params.random_seed,
                cluster_selection_algorithm=p.cluster_selection_algorithm
                ) # |> LopezAcosta2019.clean_binary_floes   # The clean-binary-floes algorithm seems pretty agressive to me.
    apply_landmask!(kmeans_result, coastal_buffer);
    nothing
end
zoom = (1000:3000, 1000:3000)
# fc_masked[zoom...]
mosaicview(Gray.(kmeans_result)[zoom...],
    Gray.(ice_water_discrim)[zoom...], nrow=1)

begin
    # The clean binary floes method has an aggressive fill_holes algorithm. Potentially merging with the
    # ice brightness threshold can prevent some of the interstitial water areas from being filled.

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = LopezAcosta2019.segmentation_B(preproc_gray, cloud_mask, kmeans_result);

    @info "tiled watershed transform on the intersect"
    @time begin
    w_merged = FSPipeline.watershed_transform(
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



## Seg B
# I'd like to update the seg B code to be more descriptive, based on what I was doing in the streamlining.
# Same with Seg F.

using IceFloeTracker
using Images
### Case tests ###
dataset = Watkins2026Dataset(; ref="v0.1")
case = first(filter(c -> (c.case_number == 11 && c.satellite == "aqua"), dataset))

segment = FSPipeline.Segment(
                cloud_mask_algorithm=LopezAcostaCloudMask(),
                coastal_buffer_structuring_element=strel_box((5,5)),
                tile_size_pixels=400,
                adapthisteq_params=(nbins=256, rblocks=8, cblocks=8, clip=1),
                )

results = segment(modis_truecolor(case), modis_falsecolor(case), modis_landmask(case))
results_old = LopezAcosta2019.Segment()(
    RGB.(modis_truecolor(case)),
    RGB.(modis_falsecolor(case)),
     RGB.(modis_landmask(case)))
begin
    cview = view_seg_random(results)
    labels = labels_map(results)
    overlay = deepcopy(modis_truecolor(case))
    overlay[labels .> 0] .= cview[labels .> 0]

    cview2 = view_seg_random(results_old)
    labels2 = labels_map(results_old)
    overlay2 = deepcopy(modis_truecolor(case))
    overlay2[labels2 .> 0] .= cview2[labels2 .> 0]
    overlay2

mosaicview(
    overlay,
    overlay2, 
    Gray.((labels .> 0) .!= (labels2 .> 0)), nrow=1)
end