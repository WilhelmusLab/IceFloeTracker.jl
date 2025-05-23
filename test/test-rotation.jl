using DataFrames
using Dates
using IceFloeTracker: get_rotation_measurements, _add_suffix

@testset "rotation" begin
    masks = Dict(
        0 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        15 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0
            0 1 1 1 1 1 1 0 0 0 1 1 1 1 1 0
            0 0 1 1 1 1 1 0 0 0 0 0 0 0 1 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        30 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 0
            0 0 0 0 0 1 0 0 0 0 0 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        45 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0
            0 0 0 0 1 1 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        60 => Bool[
            0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        75 => Bool[
            0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        90 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        105 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        120 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        135 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 1 1 1 0 0 1 0 0 0 0
            0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        150 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 1 1 1 1 0 1 1 1 0 0 0
            0 0 0 0 1 1 1 1 0 0 0 1 0 0 0 0
            0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        165 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0
            0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            0 1 1 1 1 1 0 0 0 0 1 1 1 1 0 0
            0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        180 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        195 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 1 0 0 0 0 0 0 0 1 1 1 1 1 0 0
            0 1 1 1 1 1 0 0 0 1 1 1 1 1 1 0
            0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        210 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 0 0 0 0 0 1 0 0 0 0 0
            0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        225 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 1 1 0 0 0 0
            0 0 0 0 0 0 1 1 0 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        240 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0
            0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0],
        255 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
            0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0
            0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0],
        270 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        285 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        300 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        315 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0
            0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0
            0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        330 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0
            0 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0
            0 0 0 1 1 1 0 1 1 1 1 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
            1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        345 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
            0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0
            0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
        360 => Bool[
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],
    )
    @testset "utility add suffix" begin
        df = DataFrame([
            (time=DateTime("2020-01-12T12:00:00"), mask=masks[15]),
            (time=DateTime("2020-01-13T12:00:00"), mask=masks[30])
        ])
        df_suffix1 = _add_suffix("1", df)
        @test "time1" ∈ names(df_suffix1)
        @test "mask1" ∈ names(df_suffix1)

    end

    @testset "single rows" begin

        @testset "simple case" begin

            df = DataFrame([
                (time=DateTime("2020-01-12T12:00:00"), mask=masks[0]),
                (time=DateTime("2020-01-13T12:00:00"), mask=masks[90])
            ])

            result = get_rotation_measurements(df[1, :], df[2, :];
                image_column=:mask, time_column=:time)

            @test result.theta_rad ≈ π / 2.0
            @test result.dt_sec == 86400
            @test result.row1 === df[1, :]
            @test result.row2 === df[2, :]

        end
        @testset "sampled cases" begin
            function check_single_rotation_measurements(
                time1, Δt_sec::Int64, mask1, mask2, Δθ_deg; Δθ_deg_tolerance=5.1
            )
                df = DataFrame([(time=time1, mask=mask1), (time=time1 + Second(Δt_sec), mask=mask2)])

                result = get_rotation_measurements(df[1, :], df[2, :];
                    image_column=:mask, time_column=:time)

                @test deg2rad(Δθ_deg - Δθ_deg_tolerance) <= result.theta_rad <= deg2rad(Δθ_deg + Δθ_deg_tolerance)

                @test result.dt_sec ≈ Δt_sec

                @test result.omega_rad_per_sec ≈ result.theta_rad / result.dt_sec

                @test result.row1 == df[1, :]
                @test result.row2 == df[2, :]

            end

            check_single_rotation_measurements(
                DateTime("2020-01-12T12:00:00"),
                1,
                masks[360],
                masks[0],
                0;
                Δθ_deg_tolerance=0.1,
            )
            check_single_rotation_measurements(
                DateTime("2020-01-12T12:00:00"),
                3600,
                masks[15],
                masks[0],
                -15
            )
            check_single_rotation_measurements(
                DateTime("2003-01-12T12:00:00"),
                7200,
                masks[0],
                masks[15],
                15
            )
            check_single_rotation_measurements(
                DateTime("2003-01-12T12:00:00"),
                86400,
                masks[0],
                masks[30],
                30.0
            )
            check_single_rotation_measurements(
                DateTime("2003-01-12T12:00:00"),
                864_000,
                masks[0],
                masks[30],
                30
            )
            check_single_rotation_measurements(
                DateTime("2003-01-12T12:00:00"),
                86_400 + 3600,
                masks[0],
                masks[30],
                30
            )
        end
    end

    @testset "subsets" begin
        @testset "simple case" begin

            df = DataFrame([
                (time=DateTime("2020-01-12T12:00:00"), mask=masks[0]),
                (time=DateTime("2020-01-12T13:00:00"), mask=masks[45]),
            ])

            result = DataFrame(get_rotation_measurements(df[1, :], df; image_column=:mask, time_column=:time))

            @test nrow(result) == nrow(df)
        end

        @testset "measurements match the input dataframe ordering" begin

            time = DateTime("1990-01-12T12:00:00")
            df = DataFrame([
                (time=time, mask=masks[0]),
                (time=time + Day(1), mask=masks[90]),
            ])

            result = DataFrame(get_rotation_measurements(df[2, :], df; image_column=:mask, time_column=:time))

            @test nrow(result) == nrow(df)
            @test result[1, :dt_sec] ≈ 86400.0
            @test result[1, :theta_rad] ≈ π / 2
            @test result[1, :omega_rad_per_sec] ≈ π / 2 / 24 / 3600

            @test result[2, :dt_sec] ≈ 0
            @test result[2, :theta_rad] ≈ 0
            @test isnan(result[2, :omega_rad_per_sec])  # because dt_sec == 0

            @test result[2, :row1] == result[2, :row2]
        end

    end


    @testset "full dataframe" begin
        @testset "include additional source columns" begin

            df = DataFrame([
                (id=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                (id=1, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra"),
            ])

            result = get_rotation_measurements(df; id_column=:id, image_column=:mask, time_column=:time)

            @test "satellite1" ∈ names(result)
            @test "satellite2" ∈ names(result)
            @test "mask1" ∈ names(result)
            @test "mask2" ∈ names(result)

        end
        @testset "rotation rates for an unambiguous angle" begin

            time = DateTime("2020-01-12T12:00:00")
            df = DataFrame([
                (id=1, time=time, mask=masks[0], satellite="aqua"),
                (id=1, time=time + Day(1), mask=masks[90], satellite="terra"),
            ])

            result = get_rotation_measurements(df; id_column=:id, image_column=:mask, time_column=:time)
            @test result[1, :theta_deg] ≈ 90
            @test result[1, :theta_rad] ≈ π / 2
            @test result[1, :omega_deg_per_day] ≈ 90
            @test result[1, :omega_rad_per_day] ≈ π / 2
            @test result[1, :omega_rad_per_sec] ≈ π / 2 / 24 / 3600

        end
        @testset "include additional measurement columns" begin

            df = DataFrame([
                (id=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                (id=1, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra"),
            ])

            result = get_rotation_measurements(df; id_column=:id, image_column=:mask, time_column=:time)

            @test all(result[!, :omega_rad_per_day] .≈ result[!, :omega_rad_per_sec] * 3600.0 * 24.0)
            @test all(result[!, :theta_deg] .≈ rad2deg.(result[!, :theta_rad]))
            @test all(result[!, :omega_deg_per_day] .≈ rad2deg.(result[!, :omega_rad_per_sec]) * 3600.0 * 24.0)

        end

        @testset "longer sequences" begin

            kwargs = (id_column=:id, image_column=:mask, time_column=:time)
            results = get_rotation_measurements(
                DataFrame([
                    (id=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                    (id=1, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra")]); kwargs...)
            @test nrow(results) == 1

            results = get_rotation_measurements(DataFrame([
                    (id=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                    (id=1, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra"),
                    (id=1, time=DateTime("2020-01-13T12:10:00"), mask=masks[30], satellite="aqua"),
                ]); kwargs...)
            @test nrow(results) == 3


            results = get_rotation_measurements(DataFrame([
                    (id=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                    (id=1, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra"),
                    (id=1, time=DateTime("2020-01-13T12:10:00"), mask=masks[30], satellite="aqua"),
                    (id=1, time=DateTime("2020-01-13T13:00:00"), mask=masks[30], satellite="terra"),
                ]); kwargs...)
            @test nrow(results) == 6

            results = get_rotation_measurements(DataFrame([
                    (id=1, obsid=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                    (id=1, obsid=2, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra"),
                    (id=1, obsid=3, time=DateTime("2020-01-13T12:10:00"), mask=masks[30], satellite="aqua"),
                    (id=1, obsid=4, time=DateTime("2020-01-13T13:00:00"), mask=masks[30], satellite="terra"),
                    (id=1, obsid=5, time=DateTime("2020-01-14T11:50:00"), mask=masks[45], satellite="aqua"),
                    (id=1, obsid=6, time=DateTime("2020-01-14T13:01:00"), mask=masks[90], satellite="terra"),
                ]); kwargs...)
            @test nrow(results) == 11

            @test nrow(subset(results, :obsid1 => ByRow(==(1)))) == 3  # obsid=1 is the starting position for three comparisons
            @test nrow(subset(results, :obsid2 => ByRow(==(1)))) == 0  # obsid=1 is never the "comparison" observation

            @test nrow(subset(results, :obsid1 => ByRow(==(6)))) == 0  # obsid=6 is never the "comparison" observation
            @test nrow(subset(results, :obsid2 => ByRow(==(6)))) == 3  # obsid=6 has 3 comparisons

            @test subset(results, :obsid1 => ByRow(==(5)))[1, :mask1] == masks[45]
            @test subset(results, :obsid2 => ByRow(==(6)))[1, :mask2] == masks[90]

            results = get_rotation_measurements(DataFrame([
                    (id=1, obsid=1, time=DateTime("2020-01-12T12:00:00"), mask=masks[0], satellite="aqua"),
                    (id=1, obsid=2, time=DateTime("2020-01-12T13:00:00"), mask=masks[15], satellite="terra"),
                    (id=1, obsid=3, time=DateTime("2020-01-13T12:10:00"), mask=masks[30], satellite="aqua"),
                    (id=1, obsid=4, time=DateTime("2020-01-13T13:00:00"), mask=masks[30], satellite="terra"),
                    (id=1, obsid=5, time=DateTime("2020-01-14T11:50:00"), mask=masks[45], satellite="aqua"),
                    (id=1, obsid=6, time=DateTime("2020-01-14T13:01:00"), mask=masks[90], satellite="terra"),
                    (id=2, obsid=7, time=DateTime("2020-01-12T12:00:00"), mask=masks[225], satellite="aqua"),
                    (id=2, obsid=8, time=DateTime("2020-01-12T13:00:00"), mask=masks[240], satellite="terra"),
                    (id=2, obsid=9, time=DateTime("2020-01-13T12:10:00"), mask=masks[225], satellite="aqua"),
                    (id=2, obsid=10, time=DateTime("2020-01-13T13:00:00"), mask=masks[195], satellite="terra"),
                    (id=2, obsid=11, time=DateTime("2020-01-14T11:50:00"), mask=masks[180], satellite="aqua"),
                    (id=2, obsid=12, time=DateTime("2020-01-14T13:01:00"), mask=masks[180], satellite="terra"),
                ]); kwargs...)
            @test nrow(results) == 22
            @test all(results[:, :id1] == results[:, :id2])

            @test nrow(subset(results, :obsid1 => ByRow(==(1)))) == 3  # obsid=1 is the "starting" observation three times
            @test nrow(subset(results, :obsid2 => ByRow(==(1)))) == 0  # obsid=1 is never the "comparison" observation
            @test nrow(subset(results, :obsid1 => ByRow(==(6)))) == 0  # obsid=6 is never the "starting" observation
            @test nrow(subset(results, :obsid2 => ByRow(==(6)))) == 3  # obsid=6 is the "comparison" observation three times

            @test nrow(subset(results, :obsid1 => ByRow(==(7)))) == 3  # obsid=7 is the "starting" observation three times
            @test nrow(subset(results, :obsid2 => ByRow(==(7)))) == 0  # obsid=7 is never the "comparison" observation
            @test nrow(subset(results, :obsid1 => ByRow(==(12)))) == 0  # obsid=12 is never the "starting" observation
            @test nrow(subset(results, :obsid2 => ByRow(==(12)))) == 3  # obsid=12 is the "comparison" observation three times

            @test subset(results, :obsid1 => ByRow(==(5)))[1, :mask1] == masks[45]
            @test subset(results, :obsid2 => ByRow(==(6)))[1, :mask2] == masks[90]
            @test subset(results, :obsid1 => ByRow(==(11)))[1, :mask1] == masks[180]
            @test subset(results, :obsid2 => ByRow(==(12)))[1, :mask2] == masks[180]


        end
    end
end
