# using FileIO
using Serialization
using DataFrames
using CSV



tracker_data = deserialize("notebooks/ellipses/tracker-trajectories-bug.jls")

for (i, df) in enumerate(tracker_data.props)
    CSV.write("notebooks/ellipses/example1/props-$i.csv", df)
end
