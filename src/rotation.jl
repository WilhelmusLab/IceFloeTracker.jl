add_suffix(s::String, df::DataFrame) = rename((x) -> String(x) * s, df)

function get_rotation_measurements(
    df::DataFrame; id_column, image_column, time_column
)
    results = []
    for row in eachrow(df)
        append!( # adds the 0 – n measurements from `get_rotation_measurements` to the results array
            results,
            get_rotation_measurements(
                row, df; id_column, image_column, time_column
            ),
        )
    end

    # Flatten the results into a single dataframe
    measurement_result_df = select(DataFrame(results), Not([:row1, :row2]))
    row1_df = add_suffix("1", DataFrame([r.row1 for r in results]))
    row2_df = add_suffix("2", DataFrame([r.row2 for r in results]))
    results_df = hcat(measurement_result_df, row1_df, row2_df)

    # Add some columns
    sec_per_hour = 3600.0
    sec_per_day = 86400.0
    results_df[!, "omega_rad_per_hour"] = results_df[!, "omega_rad_per_sec"] / sec_per_hour
    results_df[!, "omega_rad_per_day"] = results_df[!, "omega_rad_per_sec"] / sec_per_day

    results_df[!, "theta_deg"] = rad2deg.(results_df[!, "theta_rad"])
    results_df[!, "omega_deg_per_sec"] = rad2deg.(results_df[!, "omega_rad_per_sec"])
    results_df[!, "omega_deg_per_hour"] = results_df[!, "omega_deg_per_sec"] / sec_per_hour
    results_df[!, "omega_deg_per_day"] = results_df[!, "omega_deg_per_sec"] / sec_per_day

    return results_df

end

function get_rotation_measurements(
    measurement::DataFrameRow, df::DataFrame; id_column, image_column, time_column
)
    filtered_df = subset(
        df,
        id_column => ByRow(==(measurement[id_column])),
        time_column => ByRow((t) -> t < (measurement[time_column])), # only look at earlier images
        time_column => ByRow((t) -> Date((measurement[time_column]) - Day(1)) <= Date(t)), # only look at floes from the previous day or later
    )

    results = [
        get_rotation_measurements(
            earlier_measurement, measurement; image_column, time_column
        ) for earlier_measurement in eachrow(filtered_df)
    ]

    return results
end

function get_rotation_measurements(
    row1::DataFrameRow,
    row2::DataFrameRow;
    image_column,
    time_column,
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

