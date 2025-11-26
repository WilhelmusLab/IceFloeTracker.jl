"""
    callable_store()

Create a store and a callback function to add key-value pairs to the store.

Returns a `store::Dict` and a `callback::Function` which stores any kwargs passed to it in the `store`.

# Examples

Basic usage is to store values using the callback function
```julia-repl
julia> store, callback = callable_store()
julia> store
Dict{Any, Any}()

julia> callback(;foo="bar")  # echoes the updated store
Dict{Any, Any} with 1 entry:
  :foo => "bar"

julia> store  # values are available from the store object
Dict{Any, Any} with 1 entry:
  :foo => "bar"
```

A real-world use case is to extract data from a segmentation algorithm run:
```julia-repl
julia> intermediate_results, intermediate_results_callback = callable_store()
julia> data = first(Watkins2026Dataset(; ref="v0.1")());
julia> segments = LopezAcosta2019Tiling.Segment()(
    modis_truecolor(data),
    modis_falsecolor(data),
    modis_landmask(data);
    intermediate_results_callback,
)
Segmented Image with:
  labels map: 400×400 Matrix{Int64}
  number of labels: 12

julia> intermediate_data
Dict{Any, Any} with 16 entries:
  :binarized_tiling                       => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :icemask                                => Bool[1 1 … 1 1; 1 1 … 1 1; … ; 0 0 … 1 1; 0 0 … 1 1]
  :equalized_gray                         => [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :morphed_residue                        => [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :L0mask                                 => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :segmented                              => Segmented Image with:…
  :prelim_icemask2                        => [255 255 … 255 255; 255 255 … 255 255; … ; 255 255 … 255 255; 255 255 … 255 255]
  :equalized_gray_sharpened_reconstructed => [0 0 … 0 0; 0 0 … 0 0; … ; 255 255 … 255 255; 255 255 … 255 255]
  :gammagreen                             => [190.35 190.23 … 182.93 185.03; 191.68 190.6 … 185.04 192.08; … ; 163.87 173.33 … 108.02 108.18; 166.14 173.3 … 112.35 112.32]
  :segment_mask                           => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :ref_img_cloudmasked                    => RGB{N0f8}[RGB{N0f8}(0.0,0.0,0.0) RGB{N0f8}(0.0,0.0,0.0) … RGB{N0f8}(0.008,0.706,0.761) RGB{N0f8}(0.0,0.722,0.769); RGB{N0f8}(0.0,0.0,0.0) RGB{N0f8}(0.0,0.0,0.0) … RGB{N0f8}(0.039,0.702,0.784) RGB{N0f8}(0.075,0.784,0.859); … ; RGB{…
  :prelim_icemask                         => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :equalized_gray_reconstructed           => [0 0 … 0 0; 0 0 … 0 0; … ; 255 255 … 255 255; 255 255 … 255 255]
  :final                                  => Bool[0 0 … 0 0; 0 1 … 1 0; … ; 0 0 … 1 0; 0 0 … 0 0]
  :local_maxima_mask                      => [255 255 … 255 255; 255 255 … 255 255; … ; 255 255 … 255 255; 255 255 … 255 255]
  :labeled                                => [0 0 … 0 0; 0 1 … 1 0; … ; 0 0 … 9 0; 0 0 … 0 0]
```
"""
function callable_store()::Tuple{Dict,Function}
    store = Dict()
    function callback(; kwargs...)
        return merge!(store, Dict(kwargs))
    end
    return store, callback
end
