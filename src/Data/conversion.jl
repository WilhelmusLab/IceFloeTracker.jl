export to_landmask

to_landmask(img) = img .|> Gray .|> (x -> x .> 0.0)
