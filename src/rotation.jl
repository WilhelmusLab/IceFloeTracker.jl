using OrderedCollections: OrderedDict

function get_rotation_measurements(
    df::DataFrame; mask_column, time_column, additional_columns
)
    results = []
    for row in eachrow(df)
        append!( # adds the 0 – n measurements from `get_rotation_measurements` to the results array
            results,
            get_rotation_measurements(
                row, df; mask_column, time_column, additional_columns
            ),
        )
    end
    results_df = DataFrame(results)

    return results_df
end

function get_rotation_measurements(
    measurement::DataFrameRow, df::DataFrame; mask_column, time_column, additional_columns
)
    filtered_df = subset(
        df,
        :ID => ByRow(==(measurement[:ID])),
        time_column => ByRow((t) -> t < (measurement[time_column])), # only look at earlier images
        time_column => ByRow((t) -> Date((measurement[time_column]) - Day(1)) <= Date(t)), # only look at floes from the previous day or later
    )

    results = [
        get_rotation_measurements(
            earlier_measurement, measurement; mask_column, time_column, additional_columns
        ) for earlier_measurement in eachrow(filtered_df)
    ]

    return results
end

function get_rotation_measurements(
    row1::DataFrameRow,
    row2::DataFrameRow;
    mask_column,
    time_column,
    registration_function=register,
    # additional_columns=[],
)
    theta_rad = registration_function(row1[mask_column], row2[mask_column])
    theta_deg = rad2deg(theta_rad)

    dt = row2[time_column] - row1[time_column]
    dt_sec = dt / Dates.Second(1)
    dt_hour = dt / Dates.Hour(1)
    dt_day = dt / Dates.Day(1)

    omega_deg_per_sec = (theta_deg) / (dt_sec)
    omega_deg_per_hour = (theta_deg) / (dt_hour)
    omega_deg_per_day = (theta_deg) / (dt_day)

    omega_rad_per_sec = (theta_rad) / (dt_sec)
    omega_rad_per_hour = (theta_rad) / (dt_hour)
    omega_rad_per_day = (theta_rad) / (dt_day)

    additional_columns = setdiff(names(row1), [mask_column, time_column])

    result = OrderedDict()

    for colname in additional_columns
        result[String(colname)*"1"] = row1[colname]
        result[String(colname)*"2"] = row2[colname]
    end

    result = OrderedDict([
        "theta_deg" => theta_deg,
        "theta_rad" => theta_rad,
        String(time_column) * "1" => row1[time_column],
        String(time_column) * "2" => row2[time_column],
        "delta_time_sec" => dt_sec,
        "omega_deg_per_sec" => omega_deg_per_sec,
        "omega_deg_per_hour" => omega_deg_per_hour,
        "omega_deg_per_day" => omega_deg_per_day,
        "omega_rad_per_sec" => omega_rad_per_sec,
        "omega_rad_per_hour" => omega_rad_per_hour,
        "omega_rad_per_day" => omega_rad_per_day,
    ])

    for colname in [mask_column]
        result[String(colname)*"1"] = row1[colname]
        result[String(colname)*"2"] = row2[colname]
    end

    return result
end