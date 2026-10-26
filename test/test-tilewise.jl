@testsnippet Tilewise begin
    img = reshape(1:16, 4, 4)
    tiles = [(1:2, 1:2) (1:2, 3:4); (3:4, 1:2) (3:4, 3:4)]
end

@testitem "@tilewise stitches each tile's result into the canvas" setup = [Tilewise] begin
    @test @tilewise(tiles, fill(sum($img), size($img)), into = zeros(Int, 4, 4)) ==
        [14 14 46 46; 14 14 46 46; 22 22 54 54; 22 22 54 54]
end

@testitem "@tilewise leaves the canvas outside the tiles unchanged" setup = [Tilewise] begin
    corners = vec(tiles)[[1, 4]] # top-left and bottom-right
    @test @tilewise(corners, fill(sum($img), size($img)), into = fill(-1, 4, 4)) ==
        [14 14 -1 -1; 14 14 -1 -1; -1 -1 54 54; -1 -1 54 54]
end

@testitem "@tilewise slices only marked inputs" setup = [Tilewise] begin
    full = ones(Int, 4, 4)
    @test @tilewise(
        tiles, fill(sum($img) + length(full), size($img)), into = zeros(Int, 4, 4)
    ) == [30 30 62 62; 30 30 62 62; 38 38 70 70; 38 38 70 70]
end

@testitem "@tilewise slices several inputs" setup = [Tilewise] begin
    mask = falses(4, 4)
    mask[2, 1] = mask[3, 2] = mask[1, 4] = mask[4, 3] = true # one pixel per tile
    @test @tilewise(tiles, fill(only($img[$mask]), size($img)), into = zeros(Int, 4, 4)) ==
        [2 2 13 13; 2 2 13 13; 7 7 12 12; 7 7 12 12]
end

@testitem "@tilewise evaluates each marked input once" setup = [Tilewise] begin
    calls = Ref(0)
    counted(x) = (calls[] += 1; x)
    @tilewise tiles fill(sum($(counted(img))), 2, 2) into = zeros(Int, 4, 4)
    @test calls[] == 1
end

@testitem "@tilewise passes the tile's extent to the transformation" setup = [Tilewise] begin
    uneven = [(1:1, 1:4), (2:4, 1:4)]
    @test @tilewise(
        uneven, fill(10size($img, 1) + size($img, 2), size($img)), into = zeros(Int, 4, 4)
    ) == [14 14 14 14; 34 34 34 34; 34 34 34 34; 34 34 34 34]
end

@testitem "@tilewise keeps dimensions beyond the tile whole" setup = [Tilewise] begin
    rgb = reshape(1:48, 4, 4, 3)
    @test @tilewise(tiles, $rgb .+ sum($img), into = zeros(Int, 4, 4, 3)) ==
        rgb .+ [14 14 46 46; 14 14 46 46; 22 22 54 54; 22 22 54 54]
end

@testitem "@tilewise does not modify its inputs" setup = [Tilewise] begin
    src = collect(img)
    @tilewise tiles ($src .= 0) into = ones(Int, 4, 4)
    @test src == reshape(1:16, 4, 4)
end

@testitem "@tilewise accepts lazy and filtered tile collections" setup = [Tilewise] begin
    expected = [14 14 46 46; 14 14 46 46; 22 22 54 54; 22 22 54 54]
    for ts in (get_tiles(img; rblocks=2, cblocks=2), vec(tiles))
        @test @tilewise(ts, fill(sum($img), size($img)), into = zeros(Int, 4, 4)) == expected
    end
end

@testitem "@tilewise without a canvas collects one result per tile" setup = [Tilewise] begin
    @test @tilewise(tiles, sum($img)) == [14 46; 22 54]
    @test @tilewise(get_tiles(img; rblocks=2, cblocks=2), sum($img)) == [14 46; 22 54]
    @test @tilewise(vec(tiles), sum($img)) == [14, 22, 46, 54]
end

@testitem "@tilewise rejects a transformation with no marked input" setup = [Tilewise] begin
    @test_throws ArgumentError @macroexpand @tilewise tiles sum(img)
end

@testitem "@tilewise rejects an unknown option" setup = [Tilewise] begin
    @test_throws ArgumentError @macroexpand @tilewise tiles sum($img) onto = zeros(4, 4)
end

@testitem "tile_index_map paints each tile with its position" setup = [Tilewise] begin
    import IceFloeTracker.ImageUtils as IU
    @test IU.tile_index_map(tiles, (4, 4)) == [1 1 3 3; 1 1 3 3; 2 2 4 4; 2 2 4 4]
    @test IU.tile_index_map(get_tiles(img; rblocks=2, cblocks=2), (4, 4)) ==
        [1 1 3 3; 1 1 3 3; 2 2 4 4; 2 2 4 4]
end

@testitem "tile_index_map numbers a filtered collection by its own order" setup = [Tilewise] begin
    import IceFloeTracker.ImageUtils as IU
    corners = vec(tiles)[[1, 4]] # top-left and bottom-right
    @test IU.tile_index_map(corners, (4, 4)) == [1 1 0 0; 1 1 0 0; 0 0 2 2; 0 0 2 2]
end
