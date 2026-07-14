using IceFloeTracker
using Images

case_number = "009-" # Note: Hyphen prevents this from matching years, e.g. 201001

dataloc = "/Users/dwatkin2/Documents/research/manuscripts/greenland_floe_scale_dataset/ift_fram_strait_test_cases/data/modis/"
falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor")))
falsecolor_files = filter(f -> occursin(case_number, f), falsecolor_files)
truecolor_files = replace.(falsecolor_files, ("falsecolor" => "truecolor"))
landmask_file = joinpath(dataloc, "landmask.tiff")
landmask_img = Gray.(load(landmask_file)) .> 0;
tile_size=800
clip=1

segment = FSPipeline.Segment(;
    preprocessing_algorithm=FSPipeline.Preprocess(adapthisteq_params = (nbins=256, rblocks=8, cblocks=4, clip=clip)),
    tile_size_pixels=tile_size, # This makes the tiles close to square
    min_tile_ice_pixel_count=300
)

saveloc = "/Users/dwatkin2/Documents/research/manuscripts/greenland_floe_scale_dataset/ift_greenland_sea_dataset/data/ift_julia/calibration/"
# Vary the tile size, sharpening

for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    @info "Processing image " * tc_file
    @time begin
        segment_result = segment(tc_img, fc_img, landmask_img)
    end

    # save binary image
    save(joinpath(saveloc, "fspipeline.clip0"*string(clip)*".tile"*string(tile_size), # update this to automatically change the name with alternate settings
        replace(tc_file, ("truecolor.250m" => "segmented")), 
        ), Gray.(view_seg_random(segment_result)))
end

### Processing test case images for Aditi
using IceFloeTracker
using Images

dataloc = "/Users/dwatkin2/Documents/research/manuscripts/sea_ice_interannual_variability/UTRA2026/aditi/"
falsecolor_files = filter(f -> f != ".DS_Store", readdir(joinpath(dataloc, "falsecolor")))
truecolor_files = replace.(falsecolor_files, ("falsecolor" => "truecolor"))
landmask_file = joinpath(dataloc, "landmask.tiff")
land_mask = Gray.(load(landmask_file)) .> 0;
tile_size=1500
clip=2

segment = IceFloeTracker.FSPipeline.Segment(;
    preprocessing_algorithm=FSPipeline.Preprocess(adapthisteq_params = (nbins=256, rblocks=8, cblocks=4, clip=clip)),
    tile_size_pixels=tile_size, # This makes the tiles close to square
    min_tile_ice_pixel_count=300
)

cm = Watkins2026CloudMask(; band_7_threshold=0.15, band_2_threshold=0.34, opening_strel=strel_box((5, 5)), dilation_strel=strel_disk(5))
prelim_ice_mask = IceDetectionBrightnessPeaksMODIS134(band_1_min=0.3)

for (fc_file, tc_file) in zip(falsecolor_files, truecolor_files)
    tc_img = RGB.(load(joinpath(dataloc, "truecolor", tc_file)))
    fc_img = RGB.(load(joinpath(dataloc, "falsecolor", fc_file)))

    cloud_mask = cm(fc_img)
    
    # tbd: clear clouds that are fully contained within a coastal buffer
    apply_landmask!(cloud_mask, land_mask)

    n, m = size(tc_img)
    tile_size_pixels = tile_size
    tile_size_pixels > maximum([n, m]) && begin
        @warn "Tile size too large, defaulting to image size"
            tile_size_pixels = minimum([n, m])
    end
    
    (nr, nc) = round.(Int, size(tc_img) ./ tile_size_pixels)
    tiles = get_tiles(tc_img; rblocks=nr, cblocks=nc)

    joint_mask = land_mask .|| cloud_mask
    tc_masked = apply_landmask(tc_img, joint_mask)
    fc_masked = apply_landmask(fc_img, joint_mask)

    # First check for sufficient non-land and non-cloud pixels
    filtered_tiles = filter(
                t -> sum(.!joint_mask[t...]) > 300, tiles);

    ice_mask = prelim_ice_mask(tc_masked)

    ice_floes = segment(tc_img, fc_img, land_mask)
    floe_img = Gray.(view_seg_random(ice_floes))

    save(joinpath(dataloc, "ice_mask", replace(tc_file, "truecolor.250m.tiff"=>"ice_mask.png")), Gray.(ice_mask))
    save(joinpath(dataloc, "cloud_mask", replace(tc_file, "truecolor.250m.tiff"=>"cloud_mask.png")), Gray.(cloud_mask))
    save(joinpath(dataloc, "labeled_floes", replace(tc_file, "truecolor.250m.tiff"=>"labeled_floes.png")), floe_img)
end