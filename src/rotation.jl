using TimeZones: ZonedDateTime

_add_suffix(s::String, df::DataFrame) = rename((x) -> String(x) * s, df)

"""
Calculate the angle and rotation rate between observations in DataFrame `df`.

- `id_column` is the column with the ID of the image over several observations, e.g. the floe ID.
- `image_column` is the column with the image to compare, 
- `time_column` is the column with the timepoint of each observation,

Each row is compared to each other row in `df` which are:
  - for the same object ID,
  - strictly older,
  - not older than the previous day.

Returns a DataFrame with one row for each comparison,
with the angle `theta_rad`, time difference `dt_sec` and rotation rate `omega_rad_per_sec`,
and all the other values from `df`
with the column name suffix `1` for the first observation and `2` for the second.
"""
function get_rotation_measurements(
    df::DataFrame;
    id_column::Symbol,
    image_column::Symbol,
    time_column::Symbol,
    registration_function::Function=register,
)
    results = []
    for measurement in eachrow(df)
        filtered_df = subset(
            df,
            id_column => ByRow(==(measurement[id_column])), # only look at matching floes
            time_column => ByRow((t) -> t < (measurement[time_column])), # only look at earlier images
            time_column => ByRow((t) -> Date((measurement[time_column]) - Day(1)) <= Date(t)), # only look at floes from the previous day or later
        )
        new_results = [
            get_rotation_measurements(
                other_measurement, measurement; image_column, time_column, registration_function,
            ) for other_measurement in eachrow(filtered_df)
        ]
        push!(results, new_results)
    end
    flat_results = Iterators.flatten(results)

    # Flatten the results into a single dataframe
    measurement_result_df = select(DataFrame(flat_results), Not([:row1, :row2]))
    row1_df = _add_suffix("1", DataFrame([r.row1 for r in flat_results]))
    row2_df = _add_suffix("2", DataFrame([r.row2 for r in flat_results]))

    # Add some columns
    sec_per_day = 86400.0
    measurement_result_df[!, :omega_rad_per_day] .= measurement_result_df[!, :omega_rad_per_sec] * sec_per_day
    measurement_result_df[!, :theta_deg] .= rad2deg.(measurement_result_df[!, :theta_rad])
    omega_deg_per_sec = rad2deg.(measurement_result_df[!, :omega_rad_per_sec])
    measurement_result_df[!, :omega_deg_per_day] .= omega_deg_per_sec * sec_per_day

    results_df = hcat(measurement_result_df, row1_df, row2_df)

    return results_df

end

"""
Calculate the angle and rotation rate between an image in a DataFrameRow `measurement`,
and other images from a DataFrame `df`.

- `id_column` is the column with the ID of the image over several observations, e.g. the floe ID.
- `image_column` is the column with the image to compare, 
- `time_column` is the column with the timepoint of each observation,

Returns a vector of `NamedTuple`s with one entry for each comparison,
with the angle `theta_rad`, time difference `dt_sec` and rotation rate `omega_rad_per_sec`,
and the two input rows for each comparison `row1` and `row2`.
"""
function get_rotation_measurements(
    measurement::DataFrameRow,
    df::DataFrame;
    image_column::Symbol,
    time_column::Symbol,
    registration_function::Function=register,
)
    results = [
        get_rotation_measurements(
            other_measurement, measurement; image_column, time_column, registration_function,
        ) for other_measurement in eachrow(df)
    ]
    return results
end

"""
Calculate the angle and rotation rate between two observations in DataFrameRows `row1` and `row2`.
`image_column` and `time_column` specify which columns to use from the DataFrameRows.
Returns a NamedTuple with the angle `theta_rad`, time difference `dt_sec` and rotation rate `omega_rad_per_sec`,
and the two input rows.
"""
function get_rotation_measurements(
    row1::DataFrameRow,
    row2::DataFrameRow;
    image_column::Symbol,
    time_column::Symbol,
    registration_function::Function=register,
)
    result = get_rotation_measurements(
        row1[image_column],
        row2[image_column],
        row1[time_column],
        row2[time_column];
        registration_function
    )

    combined_result = merge(result, (; row1, row2))
    return combined_result
end

"""
Calculate the angle and rotation rate between two images `image1` and `image2` at times `time1` and `time2`.
Returns a NamedTuple with the angle `theta_rad`, time difference `dt_sec` and rotation rate `omega_rad_per_sec`.
"""
function get_rotation_measurements(
    image1::AbstractArray,
    image2::AbstractArray,
    time1::T,
    time2::T;
    registration_function::Function=register,
) where {T<:Union{ZonedDateTime,DateTime}}
    theta_rad::Float64 = registration_function(image1, image2)
    dt_sec::Float64 = (time2 - time1) / Dates.Second(1)
    omega_rad_per_sec = theta_rad / dt_sec
    result = (; theta_rad, dt_sec, omega_rad_per_sec)
    return result
end
