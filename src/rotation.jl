add_suffix(s::String, df::DataFrame) = rename((x) -> String(x) * s, df)

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
    registration_function=register,
)
    results = []
    for row in eachrow(df)
        append!( # adds the 0 – n measurements from `get_rotation_measurements` to the results array
            results,
            get_rotation_measurements(
                row, df; id_column, image_column, time_column, registration_function
            ),
        )
    end

    # Flatten the results into a single dataframe
    measurement_result_df = select(DataFrame(results), Not([:row1, :row2]))
    row1_df = add_suffix("1", DataFrame([r.row1 for r in results]))
    row2_df = add_suffix("2", DataFrame([r.row2 for r in results]))

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

`measurement` is compared to each row in the subset of `df` which are:
  - for the same object ID,
  - strictly older,
  - not older than the previous day.

Returns a vector of `NamedTuple`s with one entry for each comparison,
with the angle `theta_rad`, time difference `dt_sec` and rotation rate `omega_rad_per_sec`,
and the two input rows for each comparison `row1` and `row2`.
"""
function get_rotation_measurements(
    measurement::DataFrameRow,
    df::DataFrame;
    id_column::Symbol,
    image_column::Symbol,
    time_column::Symbol,
    registration_function=register,
)
    filtered_df = subset(
        df,
        id_column => ByRow(==(measurement[id_column])),
        time_column => ByRow((t) -> t < (measurement[time_column])), # only look at earlier images
        time_column => ByRow((t) -> Date((measurement[time_column]) - Day(1)) <= Date(t)), # only look at floes from the previous day or later
    )

    results = [
        get_rotation_measurements(
            earlier_measurement, measurement; image_column, time_column, registration_function,
        ) for earlier_measurement in eachrow(filtered_df)
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
    registration_function=register,
)
    image1::AbstractArray = row1[image_column]
    image2::AbstractArray = row2[image_column]
    theta_rad::Float64 = registration_function(image1, image2)

    datetime1 = row1[time_column]
    datetime2 = row2[time_column]
    dt = datetime2 - datetime1
    dt_sec::Float64 = dt / Dates.Second(1)

    omega_rad_per_sec = theta_rad / dt_sec

    result = (; theta_rad, dt_sec, omega_rad_per_sec, row1, row2)
    return result
end

